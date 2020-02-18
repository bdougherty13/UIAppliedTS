#' RUNS
#' @export
#' @inherit TSA::runs

RUNS<-function (x, k = 0)

{
  pruns <- function(r, n1, n2) {
    r1 <- r
    l1 <- 2
    if (n1 == n2)
      l2 <- 2 * n1
    if (n1 != n2)
      l2 <- 2 * min(n1, n2) + 1
    f <- seq(2, l2, 2)
    g1 <- seq(3, l2, 2)
    g2 <- g1
    pdf <- 0 * (1:l2)
    f[1] <- 2
    g1[1] <- n1 - 1
    g2[1] <- n2 - 1
    pdf[2] <- f[1]
    pdf[3] <- g1[1] + g2[1]
    if (l2 > 4) {
      for (i in seq(4, l2, 2)) {
        r <- (i - 2)/2
        f[r + 1] <- (n1 - r) * (n2 - r)/r/r * f[r]
        pdf[i] <- f[r + 1]
      }
    }
    if (l2 > 5) {
      for (i in seq(5, l2, 2)) {
        r <- (i - 3)/2
        g1[r + 1] <- (n1 - r - 1) * (n2 - r)/(r + 1)/r *
          g1[r]
        g2[r + 1] <- (n2 - r - 1) * (n1 - r)/(r + 1)/r *
          g2[r]
        pdf[i] <- g1[r + 1] + g2[r + 1]
      }
    }
    pdf <- pdf/sum(pdf)
    mu <- 1 + 2 * n1 * n2/(n1 + n2)
    if (r1 <= mu)
      pvalue <- sum(pdf[(1:l2) <= r1])
    if (r1 > mu)
      pvalue <- sum(pdf[(1:l2) >= r1])
    if (pvalue > 0.5)
      pvalue <- 1 - pvalue
    pvalue <- 2 * pvalue
    list(expected = mu, pvalue = signif(pvalue, 3))
  }
  x <- 1 * (x <= k)
  n1 <- sum(x)
  n2 <- length(x) - n1
  if (n1 * n2 == 0)
    return(list(pvalue = -1, expected.runs = 1 + 2 * n1 *
                  n2/(n1 + n2), n1 = n1, n2 = n2, k = k))
  r <- 1
  s <- x[1]
  for (i in 2:length(x)) {
    if (x[i] == s)
      next
    r <- r + 1
    s <- x[i]
  }
  res <- pruns(r, n1, n2)
  list(pvalue = res$pvalue, observed.runs = r, expected.runs = res$expected,
       n1 = n1, n2 = n2, k = k)
}
