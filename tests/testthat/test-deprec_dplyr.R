# MIT License
#  
#  Copyright (c) 2023 dplyr authors
#  
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#  
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#  
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
#

library(dplyr)


test_that("can select negatively (#2519)", {
  withr::local_options(lifecycle_verbosity = "quiet")

  expect_identical(monocle3::select_se(mtcars, ~ -cyl), mtcars[-2])
})  
    
test_that("select yields proper names", {
  withr::local_options(lifecycle_verbosity = "quiet")

  expect_identical(names(monocle3::select_se(mtcars, ~ cyl:hp)), c("cyl", "disp", "hp"))
})  
    

df <- tibble(
  a = c(1:3, 2:3),
  b = letters[c(1:4, 4L)]
) 


test_that("monocle3::select_se() works", {
  withr::local_options(lifecycle_verbosity = "quiet")

  expect_equal(
    monocle3::select_se(df, ~ a),
    select(df, a)
  ) 
    
  expect_equal(
    monocle3::select_se(df, ~ -a),
    select(df, -a)
  ) 
    
  expect_equal(
    monocle3::select_se(df, .dots = "a"),
    select(df, a)
  ) 
    
  expect_equal(
    monocle3::select_se(df, .dots = list(quote(-a))),
    select(df, -a)
  )
  
  expect_equal(
    monocle3::select_se(df, .dots = list(~ -a)),
    select(df, -a)
  ) 
})


