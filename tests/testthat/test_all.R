test_that("The sum constraint matrix works", {
    set.seed(100)
    r <- rnorm(100)
    mr <- matrix(r,ncol = 10, nrow = 10)
    expect_that(sum(abs(c(colSums(mr),rowSums(mr)) - genA(10,10) %*% r)), is_less_than(1e-10))
})


test_that("The log-likelihood constraint matrix works", {
    n2 <- function(x,y) pnorm(x)*pnorm(y)
    a <- seq(-2,2,by = 0.1)
    m <- outer(a,a,n2)
    BB <- genB(nrow(m),ncol(m))
    rr <- BB %*% as.vector(m)
    expect_that( length(which(rr > 0)), equals(0))
})

test_that("Recovering the matrix given the constraints works", {
    set.seed(100)
    m <- matrix(runif(15), nrow = 5)
    cr <- constraints(m)
    o <- nleqslv(as.vector(m) + runif(15)/10, slv, jacobian = jac_cr,
                 cs = cr$cs, rs = cr$rs, p = cr$p, A = cr$A, B = cr$B)
    expect_that(sum(abs(o$x - as.vector(m))), is_less_than(1e-8))
})


test_that("Trivial recovering works", {
    set.seed(101)
    m <- matrix(runif(15), nrow = 5)
    res <- recover_table(m, colSums(m), rowSums(m))
    expect_that(sum(abs(res$table-m)), is_less_than(1e-8))
})
