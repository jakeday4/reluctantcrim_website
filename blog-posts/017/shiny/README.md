# Ordinal Regression Playground — Shiny app

Interactive tool for the *Ordinal Outcomes* blog series (Part 2: *Break Things Yourself*). Embedded in the published post at `blog-posts/017/part2-break-things-yourself.html` via `<iframe>`.

## Local development

```r
shiny::runApp("app.R")
```

Required R packages: `shiny`, `ggplot2`, `dplyr`, `tidyr`, `MASS`, `ordinal`, `marginaleffects`.

## Deploy to Posit Connect Cloud

1. Sign in at [connect.posit.cloud](https://connect.posit.cloud/) with GitHub.
2. Click **Publish** → **Shiny**.
3. Pick this repo; set the working directory to `blog-posts/017/shiny`.
4. Publish. First deploy takes a few minutes while Connect Cloud installs packages.
5. Copy the public URL Connect Cloud assigns and paste it into the `app_url` line near the top of [`../part2-break-things-yourself.qmd`](../part2-break-things-yourself.qmd).

## Updating

Push changes to `app.R` on the configured branch. Connect Cloud rebuilds automatically.
