# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

css_statistic <- function(y) {
    .Call(`_variance_css_statistic`, y)
}

BS_rcpp <- function(s, e, y, penality = 1.358) {
    .Call(`_variance_BS_rcpp`, s, e, y, penality)
}

penalty_fun <- function(n, params) {
    .Call(`_variance_penalty_fun`, n, params)
}

variance <- function(v, mu) {
    .Call(`_variance_variance`, v, mu)
}

cost_function <- function(v) {
    .Call(`_variance_cost_function`, v)
}

slicing <- function(v, X, Y) {
    .Call(`_variance_slicing`, v, X, Y)
}

argmin <- function(v) {
    .Call(`_variance_argmin`, v)
}

OP_cpp <- function(data, params = 1L) {
    .Call(`_variance_OP_cpp`, data, params)
}

PELT_cpp <- function(data, params = 1L, K = 0L) {
    .Call(`_variance_PELT_cpp`, data, params, K)
}

