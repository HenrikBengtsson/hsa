## sumsq(x) is a faster version of sum(x^2)
sumsq <- function(x) {
  .Call(C_sumsq, x)
}

## sumprod(x, y) is a faster version of sum(x * y)
sumprod <- function(x, y) {
  .Call(C_sumprod, x, y)
}

## negCDbeta(C, D, beta) is a faster version of - C * D ^ beta
negCDbeta <- function(C, D, beta) {
  .Call(C_negCDbeta, C, D, beta)
}

## PERFORMANCE: dist_matrix(x) is a faster version of as.matrix(dist(x)),
## because it avoids the overhead from S3 method dispatching, handling of
## non-needed attributes etc.  Moreover, dist_matrix(x, square = TRUE) is
## avoids internal duplication of the distance matrix.
#' @importFrom stats dist
dist_matrix <- function(x, square = FALSE) {
  x <- dist(x)
  size <- attr(x, "Size")
  df <- matrix(0, nrow = size, ncol = size)
  if (square) x <- x ^ 2
  df[row(df) > col(df)] <- x
  df + t(df)
}

## PERFORMANCE: log_det(x) is a faster version of log(det(x)), because
## it avoids overhead from S3 method dispatching and skips an internal
## log(exp(t)) step.
log_det <- function(x) {
  z <- determinant.matrix(x, logarithm = TRUE)
#  d <- c(z$sign * exp(z$modulus))
#  log(d)
  if (z$sign < 0) stop("Log-determinant: NaN")
  c(z$modulus)
}

## PERFORMANCE: Avoid overhead from S3 dispatch on solve()
solve <- local({
  bs <- list()
  function(a, b = NULL) {
    if (is.null(b)) {
      n <- nrow(a)
      ## Memoization of 'b'
      if (n <= length(bs)) b <- bs[[n]]
      if (is.null(b)) {
        b <- diag(1, nrow = n)
        bs[[n]] <<- b
      }
    }
    solve.default(a, b)
  }
})

## PERFORMANCE: Avoid overhead from setting dimnames in cbind() and rbind()
cbind <- local({
  base_cbind <- base::cbind
  function(...) base_cbind(..., deparse.level = 0L)
})

rbind <- local({
  base_rbind <- base::rbind
  function(...) base_rbind(..., deparse.level = 0L)
})


## AD HOC: Trick cstruct1.R code to write files with 12 digits
## (still plenty) instead of 15 digits for easier 'diff' comparisons
#' @importFrom utils write.table
write_tsv <- function(x, ..., row.names = FALSE, col.names = FALSE,
                      sep = "\t", eol = "\n", digits = 12L) {
  for (kk in seq_along(x)) {
    if (is.double(x[[kk]])) x[[kk]] <- round(x[[kk]], digits = digits)
  }
  write.table(x, row.names = row.names, col.names = col.names,
              sep = sep, eol = eol, ...)
}