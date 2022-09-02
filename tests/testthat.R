library(testthat)
library(monocle3)

options(testthat.progress.max_fails=1000000)
Sys.setenv('TESTTHAT_MAX_FAILS' = Inf)

test_check("monocle3")

