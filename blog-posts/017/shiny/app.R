library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS)
library(ordinal)
library(marginaleffects)
library(patchwork)

# ---- Helpers ----------------------------------------------------------

gen_baseline_probs <- function(shape, K) {
  x <- 1:K
  p <- switch(
    shape,
    "Symmetric / Normal-like"    = dnorm(seq(-1.5, 1.5, length.out = K)),
    "Uniform (flat)"             = rep(1, K),
    "Bimodal (polarized)"        = {
      v <- rep(0.1, K); v[1] <- 1; v[K] <- 1; v
    },
    "Skewed toward Satisfied"    = (x / K)^2,
    "Skewed toward Dissatisfied" = ((K - x + 1) / K)^2,
    "Zero-inflated / J-shape"    = exp(-0.8 * (x - 1))
  )
  p / sum(p)
}

# Treated-arm category probabilities given baseline cell probs, latent shift β,
# effect-type, and treated-arm latent SD ratio.
#  effect_type ∈ {"constant", "asymmetric", "polarize"}
#  scale_ratio: latent SD of treated arm / latent SD of control arm. 1 = homosc.
treat_probs <- function(baseline, beta, effect_type, scale_ratio = 1) {
  K   <- length(baseline)
  tau <- qnorm(pmin(pmax(cumsum(baseline)[-K], 1e-4), 1 - 1e-4))

  if (effect_type == "polarize") {
    cum_pos <- c(pnorm((tau - beta) / scale_ratio), 1)
    cum_neg <- c(pnorm((tau + beta) / scale_ratio), 1)
    return(0.5 * diff(c(0, cum_pos)) + 0.5 * diff(c(0, cum_neg)))
  }

  beta_k <- if (effect_type == "asymmetric") {
    beta * seq(2, 0, length.out = K - 1)
  } else {
    rep(beta, K - 1)
  }

  cum_treat <- c(pnorm((tau - beta_k) / scale_ratio), 1)
  diff(c(0, cum_treat))
}

# Simulate a two-arm ordinal dataset, with optional binary subgroup covariate
# whose latent treatment effect differs from the marginal β by ±interaction_mag/2.
simulate_ab <- function(n_per, K, shape, beta, effect_type,
                        scale_ratio = 1, seed = 1138,
                        covariate = FALSE, interaction_mag = 0) {
  set.seed(seed)
  baseline <- gen_baseline_probs(shape, K)

  if (!covariate) {
    treated <- treat_probs(baseline, beta, effect_type, scale_ratio)
    ctrl <- sample(1:K, n_per, replace = TRUE, prob = baseline)
    trt  <- sample(1:K, n_per, replace = TRUE, prob = treated)
    return(tibble(
      treatment = factor(rep(c("Model A", "Model B"), each = n_per),
                         levels = c("Model A", "Model B")),
      response  = c(ctrl, trt)
    ) %>% mutate(response_ord = factor(response, levels = 1:K, ordered = TRUE)))
  }

  beta_low  <- beta + interaction_mag / 2
  beta_high <- beta - interaction_mag / 2
  n_half    <- n_per %/% 2
  n_other   <- n_per - n_half

  treated_low  <- treat_probs(baseline, beta_low,  effect_type, scale_ratio)
  treated_high <- treat_probs(baseline, beta_high, effect_type, scale_ratio)

  ctrl_low  <- sample(1:K, n_half,  replace = TRUE, prob = baseline)
  ctrl_high <- sample(1:K, n_other, replace = TRUE, prob = baseline)
  trt_low   <- sample(1:K, n_half,  replace = TRUE, prob = treated_low)
  trt_high  <- sample(1:K, n_other, replace = TRUE, prob = treated_high)

  tibble(
    treatment = factor(rep(c("Model A", "Model B"), each = n_per),
                       levels = c("Model A", "Model B")),
    subgroup  = factor(c(rep("low", n_half), rep("high", n_other),
                         rep("low", n_half), rep("high", n_other)),
                       levels = c("low", "high")),
    response  = c(ctrl_low, ctrl_high, trt_low, trt_high)
  ) %>% mutate(response_ord = factor(response, levels = 1:K, ordered = TRUE))
}

fit_models <- function(sim, fit_cs = FALSE, fit_scale = FALSE, covariate = FALSE) {
  if (covariate) {
    f_ols <- response ~ treatment * subgroup
    f_pol <- response_ord ~ treatment * subgroup
    f_cs_nominal <- ~ treatment * subgroup
  } else {
    f_ols <- response ~ treatment
    f_pol <- response_ord ~ treatment
    f_cs_nominal <- ~ treatment
  }

  m_ols <- lm(f_ols, data = sim)
  m_ord <- tryCatch(
    polr(f_pol, data = sim, method = "probit", Hess = TRUE),
    error = function(e) NULL
  )
  m_cs <- if (fit_cs) {
    tryCatch(
      clm(response_ord ~ 1, nominal = f_cs_nominal, data = sim, link = "probit"),
      error = function(e) NULL
    )
  } else NULL
  # Heteroscedastic PO probit (clm with a scale term). Gated to the
  # no-covariate, no-CS scenario per the UI conditionalPanel.
  m_scale <- if (fit_scale && !fit_cs && !covariate) {
    tryCatch(
      clm(response_ord ~ treatment, scale = ~ treatment,
          data = sim, link = "probit"),
      error = function(e) NULL
    )
  } else NULL
  list(ols = m_ols, ord = m_ord, cs = m_cs, scale = m_scale)
}

# Extract the interaction-term p-value from an lm or polr fit.
# Returns NA if there's no interaction term (e.g., no-covariate model) or
# if the model is NULL.
interaction_p <- function(model) {
  if (is.null(model)) return(NA_real_)
  cf <- summary(model)$coefficients
  rn <- rownames(cf)
  ix <- grepl(":", rn) & grepl("treatment", rn) & grepl("subgroup", rn)
  if (!any(ix)) return(NA_real_)
  if ("Pr(>|t|)" %in% colnames(cf)) return(unname(cf[ix, "Pr(>|t|)"][1]))
  if ("t value"  %in% colnames(cf)) {
    tv <- cf[ix, "t value"][1]
    return(2 * pnorm(-abs(tv)))
  }
  NA_real_
}

interaction_coef <- function(model) {
  if (is.null(model)) return(NA_real_)
  cf <- summary(model)$coefficients
  rn <- rownames(cf)
  ix <- grepl(":", rn) & grepl("treatment", rn) & grepl("subgroup", rn)
  if (!any(ix)) return(NA_real_)
  est_col <- if ("Estimate" %in% colnames(cf)) "Estimate" else "Value"
  unname(cf[ix, est_col][1])
}

predicted_probs <- function(models, K, covariate = FALSE) {
  newdat <- if (covariate) {
    expand_grid(
      treatment = factor(c("Model A", "Model B"), levels = c("Model A", "Model B")),
      subgroup  = factor(c("low", "high"), levels = c("low", "high"))
    )
  } else {
    tibble(treatment = factor(c("Model A", "Model B"), levels = c("Model A", "Model B")))
  }

  group_keys <- if (covariate) c("treatment", "subgroup") else c("treatment")

  # PO probit predictions
  if (!is.null(models$ord)) {
    raw <- predict(models$ord, newdata = newdat, type = "probs")
    if (is.matrix(raw) || is.data.frame(raw)) {
      preds_ord <- as.data.frame(raw) %>%
        setNames(as.character(1:K)) %>%
        bind_cols(newdat %>% mutate(across(everything(), as.character))) %>%
        pivot_longer(as.character(1:K), names_to = "response", values_to = "p_ord") %>%
        mutate(response = as.integer(response))
    } else {
      # vector when newdat has 1 row — shouldn't happen here, but defensive
      preds_ord <- bind_cols(
        newdat %>% mutate(across(everything(), as.character)),
        tibble(response = 1:K, p_ord = as.numeric(raw))
      )
    }
  } else {
    preds_ord <- newdat %>%
      mutate(across(everything(), as.character)) %>%
      crossing(response = 1:K) %>%
      mutate(p_ord = NA_real_)
  }

  # OLS-implied predictions: bin the OLS predictive normal at half-integer edges.
  ols_bin_probs <- function(mu, sig, K) {
    edges <- c(-Inf, seq(1.5, K - 0.5, by = 1), Inf)
    diff(pnorm(edges, mean = mu, sd = sig))
  }

  ols_mu  <- predict(models$ols, newdata = newdat)
  ols_sig <- summary(models$ols)$sigma

  preds_ols <- bind_rows(lapply(seq_len(nrow(newdat)), function(i) {
    row <- newdat[i, , drop = FALSE]
    tibble(
      !!!lapply(row, as.character),
      response = 1:K,
      p_ols    = ols_bin_probs(ols_mu[i], ols_sig, K)
    )
  }))

  # CS probit
  if (!is.null(models$cs)) {
    cs_grid <- if (covariate) {
      expand_grid(
        treatment    = factor(c("Model A", "Model B"), levels = c("Model A", "Model B")),
        subgroup     = factor(c("low", "high"), levels = c("low", "high")),
        response_ord = factor(1:K, levels = 1:K, ordered = TRUE)
      )
    } else {
      expand_grid(
        treatment    = factor(c("Model A", "Model B"), levels = c("Model A", "Model B")),
        response_ord = factor(1:K, levels = 1:K, ordered = TRUE)
      )
    }
    preds_cs <- tryCatch({
      cs_grid %>%
        mutate(p_cs = predict(models$cs, newdata = cs_grid, type = "prob")$fit,
               response = as.integer(as.character(response_ord))) %>%
        mutate(across(any_of(c("treatment", "subgroup")), as.character)) %>%
        dplyr::select(any_of(c("treatment", "subgroup", "response", "p_cs")))
    }, error = function(e) {
      newdat %>%
        mutate(across(everything(), as.character)) %>%
        crossing(response = 1:K) %>%
        mutate(p_cs = NA_real_)
    })
  } else {
    preds_cs <- newdat %>%
      mutate(across(everything(), as.character)) %>%
      crossing(response = 1:K) %>%
      mutate(p_cs = NA_real_)
  }

  # Heteroscedastic PO probit. Same per-cell prediction pattern as CS.
  # Only fitted (and only meaningful) under covariate = FALSE.
  if (!is.null(models$scale)) {
    scale_grid <- expand_grid(
      treatment    = factor(c("Model A", "Model B"), levels = c("Model A", "Model B")),
      response_ord = factor(1:K, levels = 1:K, ordered = TRUE)
    )
    preds_scale <- tryCatch({
      scale_grid %>%
        mutate(p_scale  = predict(models$scale, newdata = scale_grid, type = "prob")$fit,
               response = as.integer(as.character(response_ord))) %>%
        mutate(across(any_of(c("treatment", "subgroup")), as.character)) %>%
        dplyr::select(any_of(c("treatment", "subgroup", "response", "p_scale")))
    }, error = function(e) {
      newdat %>%
        mutate(across(everything(), as.character)) %>%
        crossing(response = 1:K) %>%
        mutate(p_scale = NA_real_)
    })
  } else {
    preds_scale <- newdat %>%
      mutate(across(everything(), as.character)) %>%
      crossing(response = 1:K) %>%
      mutate(p_scale = NA_real_)
  }

  full_join(preds_ord,   preds_ols,   by = c(group_keys, "response")) %>%
    full_join(preds_cs,    by = c(group_keys, "response")) %>%
    full_join(preds_scale, by = c(group_keys, "response"))
}

observed_probs <- function(sim, covariate = FALSE) {
  group_keys <- if (covariate) c("treatment", "subgroup") else c("treatment")
  sim %>%
    count(across(all_of(c(group_keys, "response")))) %>%
    group_by(across(all_of(group_keys))) %>%
    mutate(p_obs = n / sum(n)) %>%
    ungroup() %>%
    mutate(across(all_of(group_keys), as.character)) %>%
    dplyr::select(all_of(c(group_keys, "response", "p_obs")))
}

contrast_cis <- function(models, K, covariate = FALSE) {
  # Priority for the contrast source: category-specific probit (most flexible) >
  # heteroscedastic PO probit > standard PO probit. marginaleffects works on
  # all three; we pick whichever cumulative-probit variant the user has enabled.
  ci_src <- if (!is.null(models$cs))         models$cs
            else if (!is.null(models$scale)) models$scale
            else                              models$ord
  if (is.null(ci_src)) {
    cols <- if (covariate) {
      tibble(subgroup = rep(c("low", "high"), each = K),
             response = rep(1:K, 2))
    } else {
      tibble(response = 1:K)
    }
    return(cols %>% mutate(estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_))
  }
  me_type <- if (inherits(ci_src, "polr")) "probs" else "prob"
  tryCatch({
    out <- if (covariate) {
      avg_comparisons(ci_src, variables = "treatment", by = "subgroup", type = me_type)
    } else {
      avg_comparisons(ci_src, variables = "treatment", type = me_type)
    }
    out <- out %>%
      as_tibble() %>%
      mutate(response = suppressWarnings(as.integer(as.character(group))))
    if (covariate) {
      out %>% dplyr::select(subgroup, response, estimate, conf.low, conf.high) %>%
        mutate(subgroup = as.character(subgroup))
    } else {
      out %>% dplyr::select(response, estimate, conf.low, conf.high)
    }
  }, error = function(e) {
    cols <- if (covariate) {
      tibble(subgroup = rep(c("low", "high"), each = K),
             response = rep(1:K, 2))
    } else {
      tibble(response = 1:K)
    }
    cols %>% mutate(estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_)
  })
}

# ---- Fit statistics ---------------------------------------------------

# Log-likelihood of an OLS fit on a CATEGORICAL outcome.
# Integrates the fitted normal over half-integer bins to get a K-category
# predicted distribution, then sums log-probabilities of the observed
# responses. Puts OLS on the same likelihood scale as the cumulative
# probits so AIC / BIC comparisons are interpretable.
ols_cat_loglik <- function(model, data, K) {
  mu    <- predict(model, newdata = data)
  sig   <- summary(model)$sigma
  edges <- c(-Inf, seq(1.5, K - 0.5, by = 1), Inf)
  y     <- as.integer(data$response)
  p_lower <- pnorm(edges[y],     mean = mu, sd = sig)
  p_upper <- pnorm(edges[y + 1], mean = mu, sd = sig)
  sum(log(pmax(p_upper - p_lower, 1e-12)))
}

fit_stats <- function(models, sim, K) {
  n <- nrow(sim)
  rows <- list()

  ll_ols <- tryCatch(ols_cat_loglik(models$ols, sim, K), error = function(e) NA_real_)
  df_ols <- length(coef(models$ols)) + 1L
  rows[["Linear (OLS)"]] <- c(logLik = ll_ols, df = df_ols,
                              AIC = -2 * ll_ols + 2 * df_ols,
                              BIC = -2 * ll_ols + df_ols * log(n))

  if (!is.null(models$ord)) {
    ll <- as.numeric(logLik(models$ord))
    df <- attr(logLik(models$ord), "df")
    rows[["Cumulative probit (PO)"]] <-
      c(logLik = ll, df = df,
        AIC = -2 * ll + 2 * df,
        BIC = -2 * ll + df * log(n))
  }

  if (!is.null(models$cs)) {
    ll <- as.numeric(logLik(models$cs))
    df <- attr(logLik(models$cs), "df")
    rows[["Category-specific probit"]] <-
      c(logLik = ll, df = df,
        AIC = -2 * ll + 2 * df,
        BIC = -2 * ll + df * log(n))
  }

  if (!is.null(models$scale)) {
    ll <- as.numeric(logLik(models$scale))
    df <- attr(logLik(models$scale), "df")
    rows[["Heteroscedastic PO probit"]] <-
      c(logLik = ll, df = df,
        AIC = -2 * ll + 2 * df,
        BIC = -2 * ll + df * log(n))
  }

  do.call(rbind, rows) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Model")
}

# ---- Tour presets -----------------------------------------------------

PRESETS <- list(
  "Tour 1: Linear ≈ probit (clean)" = list(
    n_per = 500, K = 5, shape = "Symmetric / Normal-like",
    beta = 0.4, effect_type = "constant", scale_ratio = 1,
    covariate = FALSE, interaction_mag = 0,
    fit_cs = FALSE, fit_scale = FALSE,
    cf_mult = 3, seed = 1138),
  "Tour 2: Subgroup polarization (Part 1 analog)" = list(
    n_per = 1000, K = 5, shape = "Symmetric / Normal-like",
    beta = 0, effect_type = "constant", scale_ratio = 1.6,
    covariate = TRUE, interaction_mag = 1.6,
    fit_cs = TRUE, fit_scale = FALSE,
    cf_mult = 3, seed = 1138),
  "Tour 3: Skewed-tail compression" = list(
    n_per = 1000, K = 5, shape = "Skewed toward Satisfied",
    beta = 0.6, effect_type = "constant", scale_ratio = 1,
    covariate = FALSE, interaction_mag = 0,
    fit_cs = FALSE, fit_scale = FALSE,
    cf_mult = 3, seed = 1138),
  "Tour 4: Variance-only shift" = list(
    n_per = 1000, K = 5, shape = "Symmetric / Normal-like",
    beta = 0, effect_type = "constant", scale_ratio = 1.6,
    covariate = FALSE, interaction_mag = 0,
    fit_cs = FALSE, fit_scale = TRUE,
    cf_mult = 3, seed = 1138),
  "Tour 5: Counterfactual extrapolation (bimodal)" = list(
    n_per = 1000, K = 5, shape = "Bimodal (polarized)",
    beta = 0.5, effect_type = "constant", scale_ratio = 1,
    covariate = FALSE, interaction_mag = 0,
    fit_cs = FALSE, fit_scale = FALSE,
    cf_mult = 3, seed = 1138)
)

# ---- UI ---------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Ordinal Regression Playground"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      selectInput("preset", "Load preset",
                  choices = names(PRESETS),
                  selected = names(PRESETS)[1]),
      tags$hr(),

      selectInput("shape", "Baseline distribution shape",
                  choices = c("Symmetric / Normal-like",
                              "Uniform (flat)",
                              "Bimodal (polarized)",
                              "Skewed toward Satisfied",
                              "Skewed toward Dissatisfied",
                              "Zero-inflated / J-shape"),
                  selected = "Symmetric / Normal-like"),
      sliderInput("K", "Number of response categories",
                  min = 3, max = 7, value = 5, step = 1),
      sliderInput("n_per", "Sample size per arm",
                  min = 100, max = 2000, value = 500, step = 100),
      tags$hr(),

      selectInput("effect_type", "Effect type",
                  choices = c("Constant shift (PO)"           = "constant",
                              "Asymmetric (PO violation)"     = "asymmetric",
                              "Mean-preserving polarization"  = "polarize"),
                  selected = "constant"),
      sliderInput("beta", "Treatment effect on latent scale (SDs)",
                  min = -1.5, max = 1.5, value = 0.4, step = 0.05),
      sliderInput("scale_ratio", "Latent SD ratio (treated / control)",
                  min = 0.5, max = 2.0, value = 1.0, step = 0.05),
      checkboxInput("covariate",
                    "Include subgroup covariate (heterogeneous treatment)",
                    value = FALSE),
      conditionalPanel(
        "input.covariate == true",
        sliderInput("interaction_mag",
                    "Latent interaction magnitude (β_low − β_high)",
                    min = 0, max = 2, value = 0, step = 0.05),
        helpText("With β = 0 and a positive interaction magnitude, the marginal ",
                 "treatment effect is null, but the treatment helps one subgroup ",
                 "and hurts the other.")
      ),
      tags$hr(),

      checkboxInput("fit_cs",
                    "Also fit category-specific cumulative probit",
                    value = FALSE),
      conditionalPanel(
        "input.fit_cs == false && input.covariate == false",
        checkboxInput("fit_scale",
                      "Also fit heteroscedastic PO probit (clm with scale = ~ treatment)",
                      value = FALSE),
        helpText("The heteroscedastic PO probit fits a per-arm latent variance — ",
                 "the right tool when the treated arm's response distribution is ",
                 "more (or less) dispersed than the control arm's (Tour 4). ",
                 "Available only with the covariate off and category-specific probit unchecked.")
      ),
      sliderInput("cf_mult",
                  "Counterfactual treatment multiplier (Counterfactual tab)",
                  min = 1, max = 5, value = 1, step = 0.25),
      numericInput("seed", "Random seed", value = 1138, step = 1),
      helpText("Tip: pick a preset to load a canonical scenario, or build your own.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "Distributions",
          plotOutput("dist_plot", height = "550px"),
          helpText(
            "Bars = observed proportions. ",
            "Red triangles = OLS-implied category probabilities (treating the OLS ",
            "fit as a continuous normal and binning at 0.5-unit edges). Purple dots ",
            "= proportional-odds cumulative probit predicted probabilities. Green ",
            "diamonds (when 'Also fit category-specific cumulative probit' is on) = ",
            "CS probit's saturated fit. Orange squares (when 'Also fit heteroscedastic ",
            "PO probit' is on) = heteroscedastic PO probit, which fits a per-arm latent ",
            "variance — the right tool for the variance-only-shift scenario in Tour 4. ",
            "Gaps between the OLS triangles and the bars show where OLS is compressing ",
            "the tails or missing the shape."
          )
        ),
        tabPanel(
          "Coefficients & CIs",
          h4("Coefficients"),
          tableOutput("coef_table"),
          helpText(
            "OLS and the PO probit each summarize the treatment effect in a ",
            "single coefficient with a familiar p-value. The category-specific ",
            "probit cannot — its effect varies by threshold, so the contrast ",
            "table below shows per-cell estimates with 95% CIs instead. ",
            "For most applied decisions, those per-cell numbers (in proportions ",
            "of users by response category) are the deliverable; the OLS ",
            "coefficient and its p qualify a quantity — a scale-point shift in ",
            "the mean of an ordinal variable — that is commonly what decisions ",
            "ride on but, as we argue throughout the series, rarely what they ",
            "should ride on."
          ),
          conditionalPanel(
            "input.covariate == true",
            h4("Interaction tests"),
            uiOutput("interaction_box")
          ),
          h4("Model fit (AIC / BIC on a common categorical scale)"),
          tableOutput("fit_table"),
          helpText(
            "Log-likelihood, degrees of freedom, AIC, and BIC for each model fit ",
            "on the same simulated data. Lower AIC and BIC indicate better fit. ",
            "OLS's log-likelihood is computed on the categorical scale (integrating ",
            "the fitted normal over half-integer bins around each integer category), ",
            "not on the native Gaussian scale, so the three models are compared on a ",
            "common likelihood scale. Under proportional odds with a near-normal ",
            "baseline (e.g., Tour 1) the PO probit and OLS often score similarly; ",
            "under bimodal, skewed, or PO-violating settings the cumulative probit ",
            "family typically pulls ahead by a substantial margin."
          ),
          h4("Treatment contrast by category (with 95% CIs)"),
          tableOutput("contrast_table"),
          helpText(
            "The contrast is P(Y = k | Model B) − P(Y = k | Model A), by ",
            "response category (and subgroup, if the covariate is on). ",
            "Probit Δ uses whichever cumulative-probit variant is enabled: ",
            "the category-specific probit when its checkbox is on, the ",
            "heteroscedastic PO probit when that one is on, and the standard ",
            "PO probit otherwise. ",
            "⚠ flags categories where the OLS-implied contrast has the ",
            "opposite sign from the probit point estimate — those are the ",
            "cases where the choice of model not only changes the size of ",
            "the reported effect, but flips its direction."
          )
        ),
        tabPanel(
          "Stacked bars",
          plotOutput("stack_plot", height = "450px"),
          helpText(
            "100%-stacked bars show the conditional response distribution ",
            "in each arm (and subgroup, if the covariate is on). For ordinal ",
            "outcomes this is usually a cleaner summary than the mean."
          )
        ),
        tabPanel(
          "Counterfactual",
          plotOutput("cf_plot", height = "820px"),
          helpText(
            "What does each model predict at the observed treatment effect (1×) ",
            "and at the dialed counterfactual multiplier (presets default to 3×; ",
            "slide back toward 1× to collapse the comparison)? The layout depends ",
            "on the sidebar configuration. ",
            tags$br(), tags$br(),
            tags$strong("Covariate off:"), " four panels — OLS's predictive normal ",
            "and the predicted CDF on top, per-category bar comparisons for OLS ",
            "and the cumulative-probit family on the bottom, with observed proportions ",
            "in orange. The cumulative-probit panel uses the standard PO probit by ",
            "default and switches to the heteroscedastic PO probit when that checkbox ",
            "is on. Out-of-bounds ", '"<1"', " and ", '">K"', " columns are flagged ",
            "grey; cumulative probit always assigns 0 mass there, OLS does not. ",
            tags$br(), tags$br(),
            tags$strong("Covariate on:"), " the figure switches to a 2×2 bar grid ",
            "(one row per subgroup, one column per model). cf_mult scales each ",
            "subgroup's treatment effect separately — β_treatment for the reference ",
            "subgroup, β_treatment + β_interaction for the other — so the figure shows ",
            "how the heterogeneous effect extrapolates. PDF and CDF panels are dropped ",
            "in this mode; the per-subgroup bar contrasts carry the story."
          )
        )
      )
    )
  )
)

# ---- Server -----------------------------------------------------------

server <- function(input, output, session) {

  # When user picks a preset, push its values into the inputs. The selector
  # is a one-way "load this scenario" trigger — it stays on whichever preset
  # was last clicked, even if the user subsequently customizes the inputs.
  observeEvent(input$preset, {
    p <- PRESETS[[input$preset]]
    if (is.null(p)) return()
    updateSliderInput(session,  "n_per",           value = p$n_per)
    updateSliderInput(session,  "K",               value = p$K)
    updateSelectInput(session,  "shape",           selected = p$shape)
    updateSliderInput(session,  "beta",            value = p$beta)
    updateSelectInput(session,  "effect_type",     selected = p$effect_type)
    updateSliderInput(session,  "scale_ratio",     value = p$scale_ratio)
    updateCheckboxInput(session,"covariate",       value = p$covariate)
    updateSliderInput(session,  "interaction_mag", value = p$interaction_mag)
    updateCheckboxInput(session,"fit_cs",          value = p$fit_cs)
    updateCheckboxInput(session,"fit_scale",       value = p$fit_scale)
    updateSliderInput(session,  "cf_mult",         value = p$cf_mult)
    updateNumericInput(session, "seed",            value = p$seed)
  }, ignoreInit = TRUE)

  sim_dat <- reactive({
    simulate_ab(input$n_per, input$K, input$shape, input$beta,
                input$effect_type, input$scale_ratio, input$seed,
                covariate = input$covariate,
                interaction_mag = input$interaction_mag)
  })

  fits <- reactive({
    fit_models(sim_dat(),
               fit_cs    = input$fit_cs,
               fit_scale = isTRUE(input$fit_scale),
               covariate = input$covariate)
  })

  preds <- reactive({
    preds_df <- predicted_probs(fits(), input$K, covariate = input$covariate)
    obs_df   <- observed_probs(sim_dat(),  covariate = input$covariate)
    keys     <- if (input$covariate) c("treatment", "subgroup", "response")
                else                  c("treatment", "response")
    full_join(obs_df, preds_df, by = keys)
  })

  ci_df <- reactive({
    contrast_cis(fits(), input$K, covariate = input$covariate)
  })

  # Sign-flip flags: where OLS-implied Δ disagrees in sign with probit point estimate.
  sign_flip_cells <- reactive({
    dat <- preds()
    cis <- ci_df()
    keys <- if (input$covariate) c("subgroup", "response") else c("response")
    wide <- dat %>%
      pivot_wider(id_cols = all_of(keys),
                  names_from = treatment,
                  values_from = c(p_obs, p_ols)) %>%
      mutate(ols_delta = `p_ols_Model B` - `p_ols_Model A`)
    wide %>%
      left_join(cis, by = keys) %>%
      mutate(flip = !is.na(estimate) & !is.na(ols_delta) &
                    sign(estimate) * sign(ols_delta) < 0) %>%
      dplyr::select(all_of(keys), flip)
  })

  output$dist_plot <- renderPlot({
    dat <- preds()
    cs_active    <- input$fit_cs    && any(!is.na(dat$p_cs))
    scale_active <- isTRUE(input$fit_scale) && any(!is.na(dat$p_scale))

    bars_keys <- if (input$covariate) c("treatment", "subgroup", "response")
                 else                  c("treatment", "response")
    bars_dat <- dat %>%
      dplyr::select(all_of(c(bars_keys, "p_obs"))) %>%
      distinct()

    long_cols <- c("p_ord", "p_ols",
                   if (cs_active)    "p_cs",
                   if (scale_active) "p_scale")
    dat_long <- dat %>%
      pivot_longer(all_of(long_cols), names_to = "model", values_to = "p") %>%
      mutate(model = recode(model,
                            p_ord   = "Cumulative probit (PO)",
                            p_ols   = "Linear (OLS-implied)",
                            p_cs    = "Category-specific probit",
                            p_scale = "Heteroscedastic PO probit"))

    model_levels <- c("Linear (OLS-implied)", "Cumulative probit (PO)",
                      if (cs_active)    "Category-specific probit",
                      if (scale_active) "Heteroscedastic PO probit")
    dat_long$model <- factor(dat_long$model, levels = model_levels)

    color_map <- c("Linear (OLS-implied)"        = "#D3436E",
                   "Cumulative probit (PO)"      = "#341069",
                   "Category-specific probit"    = "#1B7837",
                   "Heteroscedastic PO probit"   = "#E27B0B")
    shape_map <- c("Linear (OLS-implied)"        = 17,
                   "Cumulative probit (PO)"      = 16,
                   "Category-specific probit"    = 18,
                   "Heteroscedastic PO probit"   = 15)

    p <- ggplot(mapping = aes(x = factor(response))) +
      geom_col(data = bars_dat, aes(y = p_obs),
               fill = "#FEC488", color = "#FEC488", width = 0.85) +
      geom_line(data = dat_long,
                aes(y = p, color = model, group = model),
                position = position_dodge(width = 0.5),
                alpha = 0.5, linewidth = 0.6) +
      geom_point(data = dat_long,
                 aes(y = p, color = model, shape = model),
                 size = 3, position = position_dodge(width = 0.5)) +
      scale_color_manual(values = color_map[model_levels]) +
      scale_shape_manual(values = shape_map[model_levels]) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "Response category", y = "Proportion of users",
           color = NULL, shape = NULL,
           title = "Observed vs. model-implied proportions") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", panel.grid = element_blank())

    if (input$covariate) {
      p <- p + facet_grid(subgroup ~ treatment, labeller = label_both)
    } else {
      p <- p + facet_wrap(~ treatment)
    }
    p
  })

  output$cf_plot <- renderPlot({
    models <- fits()
    sim    <- sim_dat()
    K      <- input$K
    mult   <- input$cf_mult

    validate(need(
      input$effect_type == "constant",
      "The Counterfactual view assumes a constant latent shift. Switch effect type to Constant shift (PO) to enable it."
    ))
    validate(need(
      !is.null(models$ord),
      "PO probit failed to converge. Counterfactual extrapolation requires a working PO probit fit."
    ))

    # Covariate-on branch: render the 2×2 subgroup-specific bar panel and
    # return early. The single-population PDF + CDF + bars layout below is
    # only for covariate-off scenarios.
    if (isTRUE(input$covariate)) {
      return(cf_subgroup_plot(models, sim, K, mult))
    }

    # ---- OLS extraction (1× and dialed×) ----
    newdat   <- tibble(treatment = factor(c("Model A", "Model B"),
                                          levels = c("Model A", "Model B")))
    ols_mus  <- predict(models$ols, newdata = newdat)
    ols_mu_A <- unname(ols_mus[1])
    ols_mu_B <- unname(ols_mus[2])
    ols_beta <- ols_mu_B - ols_mu_A
    ols_sig  <- summary(models$ols)$sigma
    ols_mu_cf <- ols_mu_A + mult * ols_beta

    ols_oob_lo_1x <- pnorm(0.5,         ols_mu_B,  ols_sig)
    ols_oob_hi_1x <- 1 - pnorm(K + 0.5, ols_mu_B,  ols_sig)
    ols_p_cat_1x  <- diff(pnorm(seq(0.5, K + 0.5, by = 1), ols_mu_B,  ols_sig))

    ols_oob_lo_cf <- pnorm(0.5,         ols_mu_cf, ols_sig)
    ols_oob_hi_cf <- 1 - pnorm(K + 0.5, ols_mu_cf, ols_sig)
    ols_p_cat_cf  <- diff(pnorm(seq(0.5, K + 0.5, by = 1), ols_mu_cf, ols_sig))

    # ---- PO CP / Heteroscedastic extraction (1× and dialed×) ----
    # When fit_scale is on and the heteroscedastic clm converged, we
    # substitute its predictions for PO CP's. cf_mult continues to scale
    # the latent location shift only; the variance ratio stays at whatever
    # the fit estimated (so the cf is "what if the treatment-effect mean
    # shift were Nx larger?", not "what if the variance ratio were Nx larger?").
    use_scale_model <- isTRUE(input$fit_scale) && !is.null(models$scale)
    cp_label <- if (use_scale_model) "Heteroscedastic" else "PO CP"

    if (use_scale_model) {
      # clm() with scale = ~ x stores thresholds in $alpha, location coefs in
      # $beta, and SCALE coefs in $zeta (note: this is different from polr,
      # where $zeta holds the thresholds).
      sc_thresholds <- unname(models$scale$alpha)
      sc_beta_loc   <- unname(models$scale$beta[["treatmentModel B"]])
      sc_beta_scl   <- unname(models$scale$zeta[["treatmentModel B"]])
      sigma_treated <- exp(sc_beta_scl)

      cp_p_cat_1x <- diff(c(0,
        pnorm((sc_thresholds - sc_beta_loc) / sigma_treated), 1))
      cp_p_cat_cf <- diff(c(0,
        pnorm((sc_thresholds - mult * sc_beta_loc) / sigma_treated), 1))
    } else {
      cp_beta       <- unname(coef(models$ord)[["treatmentModel B"]])
      cp_thresholds <- unname(models$ord$zeta)

      cp_p_cat_1x <- diff(c(0, pnorm(cp_thresholds - cp_beta),         1))
      cp_p_cat_cf <- diff(c(0, pnorm(cp_thresholds - mult * cp_beta),  1))
    }

    # ---- Observed proportions for Model B at 1× ----
    obs_B <- sim %>%
      dplyr::filter(treatment == "Model B") %>%
      count(response) %>%
      mutate(p_obs = n / sum(n)) %>%
      dplyr::select(response, p_obs) %>%
      right_join(tibble(response = 1:K), by = "response") %>%
      arrange(response) %>%
      mutate(p_obs = coalesce(p_obs, 0)) %>%
      pull(p_obs)

    # ---- Top-left: OLS predictive normal at dialed× ----
    bin_edges <- seq(1.5, K - 0.5, by = 1)
    x_resp <- seq(min(-0.5, ols_mu_cf - 3.5 * ols_sig),
                  max(K + 1.5, ols_mu_cf + 3.5 * ols_sig),
                  length.out = 600)
    ols_curve     <- tibble(x = x_resp, y = dnorm(x_resp, ols_mu_cf, ols_sig))
    ols_oob_left  <- ols_curve %>% dplyr::filter(x <= 0.5)
    ols_oob_right <- ols_curve %>% dplyr::filter(x >= K + 0.5)
    y_top_ols <- max(ols_curve$y) * 1.18

    p_ols_pdf <- ggplot(ols_curve, aes(x, y)) +
      geom_area(fill = "#FEC488", alpha = 0.5) +
      geom_area(data = ols_oob_left,  fill = "#888888", alpha = 0.55) +
      geom_area(data = ols_oob_right, fill = "#888888", alpha = 0.55) +
      geom_line(color = "#341069", linewidth = 0.7) +
      geom_vline(xintercept = bin_edges, linetype = "dashed",
                 color = "#341069", linewidth = 0.4) +
      geom_vline(xintercept = c(0.5, K + 0.5), linetype = "solid",
                 color = "#666666", linewidth = 0.7) +
      geom_vline(xintercept = ols_mu_cf, linetype = "solid",
                 color = "#341069", linewidth = 0.9) +
      annotate("label", x = ols_mu_cf, y = y_top_ols,
               label = sprintf("OLS predicted mean = %.2f", ols_mu_cf),
               size = 3.2, color = "#341069") +
      scale_x_continuous(breaks = seq_len(K), minor_breaks = NULL) +
      scale_y_continuous(limits = c(0, y_top_ols * 1.08),
                         expand = expansion(mult = c(0, 0))) +
      labs(x = "Response scale", y = "Density",
           title = sprintf("OLS-implied PDF at %g×", mult),
           subtitle = sprintf("Grey wings = predicted mass past the response-scale boundaries (%.0f%% total)",
                              100 * (ols_oob_lo_cf + ols_oob_hi_cf))) +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank(),
            plot.subtitle = element_text(size = 9))

    # ---- Top-right: Predicted CDF at dialed× (OLS smooth + PO CP step) ----
    cdf_curve <- tibble(x = x_resp, y = pnorm(x_resp, ols_mu_cf, ols_sig))
    cp_cdf_df <- tibble(
      x = c(min(x_resp), seq_len(K), max(x_resp)),
      y = c(0, cumsum(cp_p_cat_cf), 1)
    )

    cp_cdf_color <- if (use_scale_model) "#E27B0B" else "#341069"
    p_cdf <- ggplot() +
      geom_area(data = cdf_curve %>% dplyr::filter(x <= 0.5),
                aes(x, y), fill = "#888888", alpha = 0.30) +
      geom_area(data = cdf_curve %>% dplyr::filter(x >= K + 0.5),
                aes(x, y), fill = "#888888", alpha = 0.30) +
      geom_line(data = cdf_curve, aes(x, y, color = "OLS-implied"),
                linewidth = 0.85) +
      geom_step(data = cp_cdf_df, aes(x, y, color = cp_label),
                linewidth = 0.9, direction = "hv") +
      geom_vline(xintercept = c(0.5, K + 0.5), linetype = "solid",
                 color = "#666666", linewidth = 0.7) +
      geom_hline(yintercept = 1, linetype = "dotted",
                 color = "#999999", linewidth = 0.4) +
      scale_color_manual(values = setNames(c("#D3436E", cp_cdf_color),
                                           c("OLS-implied", cp_label)),
                         name = NULL) +
      scale_x_continuous(breaks = seq_len(K), minor_breaks = NULL,
                         limits = range(x_resp)) +
      scale_y_continuous(limits = c(-0.02, 1.05),
                         expand = expansion(mult = c(0, 0)),
                         breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      labs(x = "Response scale", y = "Cumulative probability",
           title = sprintf("Predicted CDF at %g×", mult),
           subtitle = sprintf("%s step reaches 1 at K; OLS drifts past.", cp_label)) +
      theme_minimal(base_size = 11) +
      theme(legend.position = c(0.83, 0.22),
            legend.background = element_rect(fill = "white", color = NA),
            legend.key.height = unit(0.6, "lines"),
            panel.grid = element_blank(),
            plot.subtitle = element_text(size = 9))

    # ---- Bottom: per-model bar panels (observed + 1× + dialed×) ----
    # When mult ≈ 1, the "1×" and "dialed×" series collapse onto the same
    # numbers, so we drop the duplicate series and show two bars per category
    # (observed truth + model @ 1×). When mult > 1 (or < 1), we show three.
    bar_levels  <- c("<1", as.character(seq_len(K)), paste0(">", K))
    axis_colors <- c("#666666", rep("#341069", K), "#666666")
    show_cf     <- abs(mult - 1) > 1e-6

    make_bar_panel <- function(p_1x, p_oob_lo_1x, p_oob_hi_1x,
                               p_cf, p_oob_lo_cf, p_oob_hi_cf,
                               obs_inscale, model_label,
                               color_bold, color_faint) {
      if (show_cf) {
        series_levels <- c(
          "Observed @ 1×",
          sprintf("%s @ 1×", model_label),
          sprintf("%s @ %g×", model_label, mult)
        )
        p_vec <- c(c(0, obs_inscale, 0),
                   c(p_oob_lo_1x, p_1x, p_oob_hi_1x),
                   c(p_oob_lo_cf, p_cf, p_oob_hi_cf))
        fill_values <- setNames(c("#FEC488", color_faint, color_bold),
                                series_levels)
        plot_title <- sprintf("%s: 1× vs %g× vs observed truth",
                              model_label, mult)
      } else {
        series_levels <- c(
          "Observed @ 1×",
          sprintf("%s @ 1×", model_label)
        )
        p_vec <- c(c(0, obs_inscale, 0),
                   c(p_oob_lo_1x, p_1x, p_oob_hi_1x))
        fill_values <- setNames(c("#FEC488", color_bold),
                                series_levels)
        plot_title <- sprintf("%s @ 1× vs observed truth", model_label)
      }
      n_series <- length(series_levels)
      bar_df <- tibble(
        response = factor(rep(bar_levels, n_series), levels = bar_levels),
        p        = p_vec,
        series   = factor(rep(series_levels, each = K + 2),
                          levels = series_levels)
      )
      ggplot(bar_df, aes(x = response, y = p, fill = series)) +
        annotate("rect", xmin = 0.5, xmax = 1.5,
                 ymin = 0, ymax = Inf, fill = "#888888", alpha = 0.18) +
        annotate("rect", xmin = K + 1.5, xmax = K + 2.5,
                 ymin = 0, ymax = Inf, fill = "#888888", alpha = 0.18) +
        geom_col(position = position_dodge(width = 0.85), width = 0.78) +
        scale_fill_manual(values = fill_values, name = NULL) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           expand = expansion(mult = c(0, 0.10)),
                           limits = c(0, 1)) +
        labs(x = NULL, y = "Predicted proportion",
             title = plot_title) +
        theme_minimal(base_size = 10) +
        theme(legend.position = "bottom",
              legend.key.size = unit(0.7, "lines"),
              legend.text  = element_text(size = 8),
              panel.grid   = element_blank(),
              plot.title   = element_text(size = 10),
              axis.text.x  = element_text(color = axis_colors,
                                          face  = c("bold", rep("plain", K), "bold"),
                                          size  = 9))
    }

    p_bar_ols <- make_bar_panel(ols_p_cat_1x, ols_oob_lo_1x, ols_oob_hi_1x,
                                ols_p_cat_cf, ols_oob_lo_cf, ols_oob_hi_cf,
                                obs_B, "OLS-implied",
                                "#D3436E", "#E89AAF")
    cp_colors <- if (use_scale_model) c("#E27B0B", "#F2B66E")
                 else                 c("#341069", "#8F7EAB")
    p_bar_cp  <- make_bar_panel(cp_p_cat_1x, 0, 0,
                                cp_p_cat_cf, 0, 0,
                                obs_B, cp_label,
                                cp_colors[1], cp_colors[2])

    (p_ols_pdf | p_cdf) / (p_bar_ols | p_bar_cp) +
      plot_layout(heights = c(1, 1.05))
  })

  # ---- Subgroup-specific counterfactual (covariate = TRUE) -----------
  # When a subgroup covariate is on and the effect type is constant, the
  # Counterfactual tab renders a 2 × 2 bar panel: one row per subgroup,
  # one column per model (OLS vs PO CP). No PDF/CDF panels — the bar
  # comparison is the point. cf_mult scales each subgroup's treatment
  # effect separately (β_treatment for the reference subgroup, β_treatment
  # + β_interaction for the other), so the counterfactual respects the
  # heterogeneous treatment structure.
  cf_subgroup_plot <- function(models, sim, K, mult) {
    ols_coefs <- coef(models$ols)
    ols_sig   <- summary(models$ols)$sigma
    ols_b_t   <- ols_coefs[["treatmentModel B"]]
    ols_b_s   <- ols_coefs[["subgrouphigh"]]
    ols_b_x   <- ols_coefs[["treatmentModel B:subgrouphigh"]]
    ols_int   <- ols_coefs[["(Intercept)"]]

    # OLS predicted means per (treatment, subgroup) cell at 1× and mult×.
    # Counterfactual scales only the treatment-related effects.
    ols_mu_A_low  <- ols_int
    ols_mu_A_high <- ols_int + ols_b_s
    ols_mu_B_low_1x  <- ols_int + ols_b_t
    ols_mu_B_high_1x <- ols_int + ols_b_s + ols_b_t + ols_b_x
    ols_mu_B_low_cf  <- ols_int + mult * ols_b_t
    ols_mu_B_high_cf <- ols_int + ols_b_s + mult * (ols_b_t + ols_b_x)

    ols_cell <- function(mu) {
      edges <- c(-Inf, seq(0.5, K + 0.5, by = 1), Inf)
      pp    <- diff(pnorm(edges, mean = mu, sd = ols_sig))
      list(p_oob_lo = pp[1],
           p_cat    = pp[2:(K + 1)],
           p_oob_hi = pp[K + 2])
    }

    # PO probit per-subgroup latent shifts.
    cp_coefs <- coef(models$ord)
    cp_thresh <- models$ord$zeta
    cp_b_t  <- cp_coefs[["treatmentModel B"]]
    cp_b_s  <- cp_coefs[["subgrouphigh"]]
    cp_b_x  <- cp_coefs[["treatmentModel B:subgrouphigh"]]

    cp_cells <- function(treatment_eff, subgroup_eff) {
      # P(Y <= k | cell) = Φ(τ_k - (treatment_eff + subgroup_eff))
      diff(c(0, pnorm(cp_thresh - (treatment_eff + subgroup_eff)), 1))
    }

    # Observed Model B proportions per subgroup
    obs_per_subgroup <- function(sg) {
      sim %>%
        dplyr::filter(treatment == "Model B", subgroup == sg) %>%
        count(response) %>%
        right_join(tibble(response = 1:K), by = "response") %>%
        arrange(response) %>%
        mutate(p = coalesce(n, 0L) / sum(coalesce(n, 0L))) %>%
        pull(p)
    }
    obs_low  <- obs_per_subgroup("low")
    obs_high <- obs_per_subgroup("high")

    bar_levels  <- c("<1", as.character(seq_len(K)), paste0(">", K))
    axis_colors <- c("#666666", rep("#341069", K), "#666666")
    show_cf     <- abs(mult - 1) > 1e-6

    make_subgroup_panel <- function(p_1x, p_oob_lo_1x, p_oob_hi_1x,
                                    p_cf, p_oob_lo_cf, p_oob_hi_cf,
                                    obs_inscale, model_label, subgroup_label,
                                    color_bold, color_faint) {
      if (show_cf) {
        series_levels <- c(
          "Observed @ 1×",
          sprintf("%s @ 1×", model_label),
          sprintf("%s @ %g×", model_label, mult)
        )
        p_vec <- c(c(0, obs_inscale, 0),
                   c(p_oob_lo_1x, p_1x, p_oob_hi_1x),
                   c(p_oob_lo_cf, p_cf, p_oob_hi_cf))
        fill_values <- setNames(c("#FEC488", color_faint, color_bold),
                                series_levels)
      } else {
        series_levels <- c(
          "Observed @ 1×",
          sprintf("%s @ 1×", model_label)
        )
        p_vec <- c(c(0, obs_inscale, 0),
                   c(p_oob_lo_1x, p_1x, p_oob_hi_1x))
        fill_values <- setNames(c("#FEC488", color_bold), series_levels)
      }
      n_series <- length(series_levels)
      bar_df <- tibble(
        response = factor(rep(bar_levels, n_series), levels = bar_levels),
        p        = p_vec,
        series   = factor(rep(series_levels, each = K + 2),
                          levels = series_levels)
      )
      ggplot(bar_df, aes(x = response, y = p, fill = series)) +
        annotate("rect", xmin = 0.5, xmax = 1.5,
                 ymin = 0, ymax = Inf, fill = "#888888", alpha = 0.18) +
        annotate("rect", xmin = K + 1.5, xmax = K + 2.5,
                 ymin = 0, ymax = Inf, fill = "#888888", alpha = 0.18) +
        geom_col(position = position_dodge(width = 0.85), width = 0.78) +
        scale_fill_manual(values = fill_values, name = NULL) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           expand = expansion(mult = c(0, 0.10)),
                           limits = c(0, 1)) +
        labs(x = NULL, y = "Predicted proportion",
             title = sprintf("%s — %s subgroup", model_label, subgroup_label)) +
        theme_minimal(base_size = 10) +
        theme(legend.position = "bottom",
              legend.key.size = unit(0.7, "lines"),
              legend.text  = element_text(size = 8),
              panel.grid   = element_blank(),
              plot.title   = element_text(size = 10),
              axis.text.x  = element_text(color = axis_colors,
                                          face  = c("bold", rep("plain", K), "bold"),
                                          size  = 9))
    }

    # Build the 4 panels
    ols_B_low_1x   <- ols_cell(ols_mu_B_low_1x)
    ols_B_low_cf   <- ols_cell(ols_mu_B_low_cf)
    ols_B_high_1x  <- ols_cell(ols_mu_B_high_1x)
    ols_B_high_cf  <- ols_cell(ols_mu_B_high_cf)

    cp_B_low_1x  <- cp_cells(cp_b_t,                 0)
    cp_B_low_cf  <- cp_cells(mult * cp_b_t,          0)
    cp_B_high_1x <- cp_cells(cp_b_t + cp_b_x,        cp_b_s)
    cp_B_high_cf <- cp_cells(mult * (cp_b_t + cp_b_x), cp_b_s)

    p_ols_low <- make_subgroup_panel(
      ols_B_low_1x$p_cat, ols_B_low_1x$p_oob_lo, ols_B_low_1x$p_oob_hi,
      ols_B_low_cf$p_cat, ols_B_low_cf$p_oob_lo, ols_B_low_cf$p_oob_hi,
      obs_low, "OLS-implied", "low", "#D3436E", "#E89AAF")
    p_ols_high <- make_subgroup_panel(
      ols_B_high_1x$p_cat, ols_B_high_1x$p_oob_lo, ols_B_high_1x$p_oob_hi,
      ols_B_high_cf$p_cat, ols_B_high_cf$p_oob_lo, ols_B_high_cf$p_oob_hi,
      obs_high, "OLS-implied", "high", "#D3436E", "#E89AAF")
    p_cp_low <- make_subgroup_panel(
      cp_B_low_1x, 0, 0, cp_B_low_cf, 0, 0,
      obs_low, "PO CP", "low", "#341069", "#8F7EAB")
    p_cp_high <- make_subgroup_panel(
      cp_B_high_1x, 0, 0, cp_B_high_cf, 0, 0,
      obs_high, "PO CP", "high", "#341069", "#8F7EAB")

    (p_ols_low | p_cp_low) / (p_ols_high | p_cp_high) +
      plot_layout(heights = c(1, 1))
  }

  output$stack_plot <- renderPlot({
    sim <- sim_dat()
    p <- ggplot(sim, aes(x = treatment, fill = factor(response))) +
      geom_bar(position = "fill") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_fill_viridis_d(option = "magma", begin = 0.15, end = 0.85,
                           name = "Response", direction = -1) +
      labs(x = NULL, y = "Proportion of users") +
      coord_flip() +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank())
    if (input$covariate) p <- p + facet_wrap(~ subgroup, labeller = label_both)
    p
  })

  output$coef_table <- renderTable({
    models <- fits()

    fmt_p <- function(p) {
      if (is.na(p)) return("—")
      if (p < 0.001) return("< 0.001")
      sprintf("%.3f", p)
    }

    # For interaction (covariate) models, the "treatmentModel B" coefficient is
    # the effect at the reference subgroup, which is misleading to display by
    # itself — see the *Interaction tests* panel for the relevant inference.
    if (input$covariate) {
      tibble(
        Model          = c("Linear (OLS)", "Cumulative probit (PO)",
                           "Category-specific probit"),
        Coefficient    = rep("(see interaction tests below)", 3),
        `p-value`      = rep("—", 3),
        Interpretation = c("scale-point shift in mean response",
                           "latent-SD shift",
                           "per-threshold shifts")
      )
    } else {
      ols_summary <- summary(models$ols)$coefficients
      ols_coef <- sprintf("%.3f", ols_summary["treatmentModel B", "Estimate"])
      ols_p    <- fmt_p(ols_summary["treatmentModel B", "Pr(>|t|)"])

      if (!is.null(models$ord)) {
        ord_summary <- summary(models$ord)$coefficients
        ord_coef <- sprintf("%.3f", ord_summary["treatmentModel B", "Value"])
        ord_p    <- fmt_p(2 * pnorm(-abs(ord_summary["treatmentModel B", "t value"])))
      } else {
        ord_coef <- "— (fit failed)"
        ord_p    <- "—"
      }

      cs_note <- if (!is.null(models$cs)) {
        "(per-threshold; see contrast table)"
      } else if (input$fit_cs) "— (fit failed)" else "—"
      cs_p <- "—"

      # Heteroscedastic PO: report location β (latent SDs) and the scale
      # ratio (treated SD / control SD) implied by the fitted scale coef.
      # clm stores location coefs in $beta and scale coefs in $zeta when
      # a scale formula is provided (different convention from polr).
      if (!is.null(models$scale)) {
        loc_b <- unname(models$scale$beta[["treatmentModel B"]])
        sc_b  <- tryCatch(
          unname(models$scale$zeta[["treatmentModel B"]]),
          error = function(e) NA_real_)
        scale_coef <- if (!is.na(sc_b)) {
          sprintf("β = %.3f, σ-ratio = %.3f", loc_b, exp(sc_b))
        } else {
          sprintf("β = %.3f (scale extraction failed)", loc_b)
        }
        scale_p <- "(see fit table)"
      } else {
        scale_coef <- if (isTRUE(input$fit_scale)) "— (fit failed)" else "—"
        scale_p    <- "—"
      }

      tibble(
        Model          = c("Linear (OLS)", "Cumulative probit (PO)",
                           "Category-specific probit",
                           "Heteroscedastic PO probit"),
        Coefficient    = c(ols_coef, ord_coef, cs_note, scale_coef),
        `p-value`      = c(ols_p, ord_p, cs_p, scale_p),
        Interpretation = c("scale-point shift in mean response",
                           "latent-SD shift",
                           "per-threshold shifts",
                           "latent-SD location shift + latent-SD scale ratio")
      )
    }
  })

  output$fit_table <- renderTable({
    fit_stats(fits(), sim_dat(), input$K) %>%
      dplyr::mutate(
        logLik = sprintf("%.1f", logLik),
        df     = as.integer(df),
        AIC    = sprintf("%.1f", AIC),
        BIC    = sprintf("%.1f", BIC)
      )
  })

  output$interaction_box <- renderUI({
    if (!input$covariate) return(NULL)
    models <- fits()
    ols_coef <- interaction_coef(models$ols)
    ols_p    <- interaction_p(models$ols)
    pol_coef <- interaction_coef(models$ord)
    pol_p    <- interaction_p(models$ord)

    fmt_p    <- function(p) if (is.na(p)) "—" else if (p < 0.001) "< 0.001" else sprintf("%.3f", p)
    fmt_coef <- function(b) if (is.na(b)) "—" else sprintf("%.3f", b)

    tags$div(
      tags$p(
        tags$strong("Linear (OLS):"),
        sprintf(" treatment × subgroup interaction = %s, p = %s.",
                fmt_coef(ols_coef), fmt_p(ols_p))
      ),
      tags$p(
        tags$strong("Cumulative probit (PO):"),
        sprintf(" treatment × subgroup interaction = %s, p = %s.",
                fmt_coef(pol_coef), fmt_p(pol_p))
      ),
      tags$p(
        tags$em(
          "Two models, two verdicts on the same data. The probit reports the ",
          "latent-scale interaction; OLS reports the response-scale interaction. ",
          "When baseline shapes are non-symmetric or near boundaries, these can ",
          "diverge sharply — even on whether an interaction exists at all."
        )
      )
    )
  })

  output$contrast_table <- renderTable({
    models <- fits()
    validate(need(
      !is.null(models$ord) || !is.null(models$cs),
      "Probit models failed to converge for these settings. Try a larger sample size, a less extreme baseline shape, or a smaller effect."
    ))

    dat   <- preds()
    cis   <- ci_df()
    flips <- sign_flip_cells()

    pct <- function(x) ifelse(is.na(x), "—", sprintf("%.0f%%", 100 * x))
    with_ci <- function(est, lo, hi) {
      if (is.na(est)) return("—")
      sprintf("%.0f%% [%.0f%%, %.0f%%]", 100 * est, 100 * lo, 100 * hi)
    }

    keys <- if (input$covariate) c("subgroup", "response") else c("response")

    contrasts_pt <- dat %>%
      pivot_wider(id_cols = all_of(keys),
                  names_from = treatment,
                  values_from = c(p_obs, p_ols)) %>%
      mutate(
        obs_delta = `p_obs_Model B` - `p_obs_Model A`,
        ols_delta = `p_ols_Model B` - `p_ols_Model A`
      ) %>%
      dplyr::select(all_of(keys), obs_delta, ols_delta)

    out <- contrasts_pt %>%
      left_join(cis,   by = keys) %>%
      left_join(flips, by = keys) %>%
      arrange(across(all_of(keys))) %>%
      mutate(
        Category            = as.character(response),
        `Observed Δ`        = pct(obs_delta),
        `OLS-implied Δ`     = ifelse(!is.na(flip) & flip,
                                     paste0(pct(ols_delta), " ⚠"),
                                     pct(ols_delta)),
        `Probit Δ [95% CI]` = mapply(with_ci, estimate, conf.low, conf.high)
      )

    if (input$covariate) {
      out %>% dplyr::select(Subgroup = subgroup, Category,
                            `Observed Δ`, `OLS-implied Δ`, `Probit Δ [95% CI]`)
    } else {
      out %>% dplyr::select(Category, `Observed Δ`, `OLS-implied Δ`,
                            `Probit Δ [95% CI]`)
    }
  })
}

shinyApp(ui, server)
