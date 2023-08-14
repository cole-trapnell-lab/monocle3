library(testthat)
library(monocle3)

# Prevent flashing text on some Linux platforms.
options(testthat..use_colours = FALSE)
# Run to end of tests.
options(testthat.progress.max_fails=1000000)
Sys.setenv('TESTTHAT_MAX_FAILS' = Inf)

test_check("monocle3")

