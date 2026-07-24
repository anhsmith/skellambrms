# Precompile the vignette.
#
# vignettes/pairedcountbrms.Rmd fits a real brms model. Building that on every
# R CMD check and every CI run would need a full Stan toolchain plus minutes of
# compilation, so it is precompiled instead: the .Rmd.orig source is knitted
# HERE, on a machine with Stan, and the resulting .Rmd carries the output as
# static text. Downstream builds only render markdown.
#
# Run this by hand after changing pairedcountbrms.Rmd.orig, then commit BOTH
# files. The knitr::knit() call executes every chunk, so the numbers in the
# committed .Rmd are real output, not transcribed.
#
#   Rscript vignettes/precompile.R
#
# The vignette produces no figures, so there is no fig.path to rewrite. If you
# add a plotting chunk, set fig.path to "figure/" and make sure the generated
# images are committed alongside the .Rmd.

stopifnot(
  "run from the package root" = file.exists("DESCRIPTION"),
  "needs a Stan backend" =
    requireNamespace("cmdstanr", quietly = TRUE) ||
    requireNamespace("rstan", quietly = TRUE)
)

message("Knitting pairedcountbrms.Rmd.orig -- this fits a model, allow a few minutes.")

old <- setwd("vignettes")
on.exit(setwd(old), add = TRUE)

knitr::knit(
  input  = "pairedcountbrms.Rmd.orig",
  output = "pairedcountbrms.Rmd"
)

message("Done. Commit both pairedcountbrms.Rmd.orig and pairedcountbrms.Rmd.")
