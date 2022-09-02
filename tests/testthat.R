library(testthat)
library(monocle3)

option(testthat.progress.max_fails=1000000)

test_check("monocle3")

