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
#
#  Deprecated SE versions of main verbs.
#
#  dplyr used to offer twin versions of each verb suffixed with an
#  underscore. These versions had standard evaluation (SE) semantics:
#  rather than taking arguments by code, like NSE verbs, they took
#  arguments by value. Their purpose was to make it possible to
#  program with dplyr. However, dplyr now uses tidy evaluation
#  semantics. NSE verbs still capture their arguments, but you can now
#  unquote parts of these arguments. This offers full programmability
#  with NSE verbs. Thus, the underscored versions are now superfluous.
#
#  Unquoting triggers immediate evaluation of its operand and inlines
#  the result within the captured expression. This result can be a
#  value or an expression to be evaluated later with the rest of the
#  argument. See `vignette("programming")` for more information.
#


lazy_deprec <- function(fun,
                        hint = TRUE,
                        env = caller_env(),
                        user_env = caller_env(2)) {
  lifecycle::deprecate_warn(
    env = env,
    user_env = user_env,
    always = TRUE
  )
}


warn_underscored_se <- function() {
  return(NULL)
  warn(paste(
    "The underscored versions are deprecated in favour of",
    "tidy evaluation idioms. Please see the documentation",
    "for `quo()` in rlang"
  ))
}


warn_text_se <- function() {
  return(NULL)
  warn("Text parsing is deprecated, please supply an expression or formula")
}


compat_lazy_se <- function(lazy, env = caller_env(), warn = TRUE) {
  # Note: warn_underscored is disabled above.
  if (warn) monocle3:::warn_underscored_se()

  if (missing(lazy)) {
    return(quo())
  }
  if (is_quosure(lazy)) {
    return(lazy)
  }
  if (is_formula(lazy)) {
    return(as_quosure(lazy, env))
  }

  out <- switch(typeof(lazy),
    symbol = ,
    language = new_quosure(lazy, env),
    character = {
      if (warn) monocle3:::warn_text_se()
      parse_quo(lazy[[1]], env)
    },
    logical = ,
    integer = ,
    double = {
      if (length(lazy) > 1) {
        warn("Truncating vector to length 1")
        lazy <- lazy[[1]]
      }
      new_quosure(lazy, env)
    },
    list =
      if (inherits(lazy, "lazy")) {
        lazy = new_quosure(lazy$expr, lazy$env)
      }
  )

  if (is_null(out)) {
    abort(sprintf("Can't convert a %s to a quosure", typeof(lazy)))
  } else {
    out
  }
}


compat_lazy_dots_se <- function(dots, env, ..., .named = FALSE) {
  if (missing(dots)) {
    dots <- list()
  }
  if (inherits(dots, c("lazy", "formula"))) {
    dots <- list(dots)
  } else {
    dots <- unclass(dots)
  }
  dots <- c(dots, list(...))

  warn <- TRUE
  for (i in seq_along(dots)) {
    dots[[i]] <- monocle3:::compat_lazy_se(dots[[i]], env, warn)
    warn <- FALSE
  }

  named <- have_name(dots)
  if (.named && any(!named)) {
    nms <- vapply(dots[!named], function(x) expr_text(get_expr(x)), character(1))
    names(dots)[!named] <- nms
  }

  names(dots) <- names2(dots)
  dots
}


# Generic select_se.
#' @export
select_se <- function(.data, ..., .dots = list()) {
  # lazy_deprec("select", hint = FALSE) # Disable warning message.
  UseMethod("select_se")
}


# select_se for data.frames.
#' @export
select_se.data.frame <- function(.data, ..., .dots = list()) {
  dots <- monocle3:::compat_lazy_dots_se(.dots, caller_env(), ...)
  dplyr::select(.data, !!!dots)
}


# select_se for grouped data frames.
#' @export
select_se.grouped_df <- function(.data, ..., .dots = list()) {
  dots <- monocle3:::compat_lazy_dots_se(.dots, caller_env(), ...)
  dplyr::select(.data, !!!dots)
}


