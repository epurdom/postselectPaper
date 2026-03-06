# postselect

R package for post-selection inference (customize this description).

## Installation

```r
# From GitHub (after you push):
# remotes::install_github("your-org/postselectPaper")
```

## Development

- **R code** goes in `R/`. Use roxygen2 comments to document functions; then run `devtools::document()` (or Ctrl+Shift+D in RStudio) to update `man/` and `NAMESPACE`.
- **Tests** go in `tests/testthat/` using testthat. Run with `devtools::test()`.
- **Data** for examples/vignettes: put raw data in `data-raw/` and document how to create `data/*.rda`; or put small data in `inst/extdata/`.
- **Vignettes** (optional): add `vignettes/` and use `devtools::build_vignettes()`.

After moving files in, run:

```r
devtools::document()
devtools::check()
```

## License

MIT
