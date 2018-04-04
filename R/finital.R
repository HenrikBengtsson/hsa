#' @importFrom stats optim rnorm
finital <- function(pbin, A0, b0, invSigma0, beta1, covmat0, mat, floglike, fdloglike) {
  N <- length(pbin)
  if (N <= 500) {
    P <- Sis(5, pbin, A0, b0, invSigma0, beta1, covmat0, mat, c(0, 0, 0), function(x, ...) -loglikelihood(x, ...))
  } else {
    if (N > 500 && N <= 2000) {
      m <- 100
    } else {
      m <- 200
    }
    index <- cbind(c(1, seq(from = m, to = N - m, by = m)), c(seq(from = m, to = N - m, by = m), N))
    index[-1, 1] <- index[-1, 1] + 1L
    lP <- lapply(1:dim(index)[1], FUN = function(x) {
      idxs <- index[x, 1]:index[x, 2]
      subinitial(pbin[idxs], A0, b0, invSigma0, beta1, covmat0[idxs, idxs, ], mat[idxs, c(1, idxs + 1L), ], floglike, fdloglike)
    })
    P <- matrix(0, nrow = N, ncol = 3L)
    P[index[1, 1]:index[1, 2], ] <- lP[[1]]
    for (i in 2:dim(index)[1]) {
      index_ri_c1 <- index[i, 1]
      index_ri_c2 <- index[i, 2]
      idxs <- 1:index_ri_c2
      P_i <- P[1:(index[i - 1, 2]), ]
      pbin_i <- pbin[idxs]
      lP_i <- lP[[i]]
      covmat0_i <- covmat0[idxs, idxs, ]
      mat_i <- mat[idxs, 1:(index_ri_c2 + 1L), ]
      theta <- optim(as.vector(rbind(diag(3), P[index[i - 1, 2], ] - P[index_ri_c1, ] + rnorm(3L, sd = 1/100))), fn = function(x) {
        piece(x, P_i, lP_i, pbin_i, A0, b0, invSigma0, beta1, covmat0_i, mat_i, floglike)
      })
      theta <- matrix(theta[["par"]], nrow = 4L, ncol = 3L)
      P[index_ri_c1:index_ri_c2, ] <- lP_i %*% theta[-4, ] + rep(1, times = index_ri_c2 - index_ri_c1 + 1L) %*% t(theta[4, ])
    }
  }
  avsmth(pbin, P = P)
}
