test_that("The sum constraint matrix works", {
    set.seed(100)
    r <- rnorm(100)
    mr <- matrix(r,ncol = 10, nrow = 10)
    expect_that(sum(abs(c(colSums(mr),rowSums(mr)) - genA(10,10) %*% r)), is_less_than(1e-10))
})


test_that("The log ratio constraint matrix works", {
    n2 <- function(x,y) pnorm(x)*pnorm(y)
    a <- seq(-2,2,by = 0.1)
    m <- outer(a,a,n2)
    BB <- genB(nrow(m),ncol(m))
    rr <- BB %*% as.vector(m)
    expect_that( length(which(rr > 0)), equals(0))
})

test_that("The inverse of fixed constraints works", {
    set.seed(12)
    m <- matrix(runif(15),nrow = 5)
    cr <- constraints(m, ratio = "fixed")
    mm <- invp1(cr$p, m[, 1], m[1, ])
    expect_that(sum(abs(mm - m)), is_less_than(1e-6))

})

test_that("The inverse of sequential constraints works", {
    set.seed(112)
    m <- matrix(runif(15),nrow = 5)
    cr <- constraints(m, ratio = "sequential")
    mm <- invp(cr$p, m[, 1], m[1, ])
    expect_that(sum(abs(mm - m)), is_less_than(1e-6))

})

test_that("Getting rows from columns and log ratio works", {
    set.seed(113)
    m <- matrix(runif(28),nrow = 7)
    cr <- constraints(m)
    row1 <- get_r(m[,1], cr$cs, cr$rs, cr$p)
    expect_that(sum(abs(row1 - m[1,-1])), is_less_than(1e-6))

})

test_that("Getting columns from rows and log ratio works", {
    set.seed(114)
    m <- matrix(runif(63),nrow = 9)
    cr <- constraints(m)
    col1 <- get_c(m[1,-1], cr$cs, cr$rs, cr$p)
    expect_that(sum(abs(col1 - m[,1])), is_less_than(1e-6))
})

test_that("The log ratio constraint with fixed first row and column  works", {
    set.seed(13)
    m <- matrix(runif(15),nrow = 5)
    BB <- genB(nrow(m),ncol(m))
    BB1 <- genB1(nrow(m),ncol(m))
    rr <- matrix(BB %*% as.vector(m),nrow = nrow(m) - 1)
    rr1 <- matrix(BB1 %*% as.vector(m),nrow = nrow(m) - 1)
    expect_that( sum(abs(rr1 - cumsum2(rr))), is_less_than(1e-6))
})

test_that("The jacobian works", {
    set.seed(115)
    m <- matrix(runif(63),nrow = 9)
    cr <- constraints(m)
    jc_t <- jac_cr(cr$ix, p = cr$p, A = cr$A)
    jc_n <- jacobian(slv, cr$ix, cs = cr$cs, rs =cr$rs, p = cr$p)
    expect_that( sum(abs(jc_t - jc_n)), is_less_than(1e-6))
})

test_that("Recovering the matrix given the constraints works", {
    set.seed(100)
    m <- matrix(runif(15), nrow = 5)
    cr <- constraints(m)
    o <- nleqslv(cr$ix + runif(length(cr$ix))/14, slv, jacobian = jac_cr,
                 cs = cr$cs, rs = cr$rs, p = cr$p, A = cr$A)

    res <- invp1(cr$p, o$x[1:5], c(o$x[1], o$x[6:7]))
    expect_that(sum(abs(res - m)), is_less_than(1e-6))
})


test_that("Trivial rescaling works", {
    set.seed(101)
    m <- matrix(runif(15), nrow = 5)
    res <- rescale_table(m, colSums(m), rowSums(m))
    expect_that(sum(abs(res$table - m)), is_less_than(1e-8))
})

test_that("Trivial recovering works with fixed ratios", {
    set.seed(102)
    mm <- matrix(runif(15), nrow = 5)
    cr <-  constraints(mm, ratio = "fixed")
    res <- recover_table(cr$p, colSums(mm), rowSums(mm), ratio = "fixed")
    expect_that(sum(abs(res$table - mm)), is_less_than(1e-8))
})


test_that("Trivial recovering works with sequential ratios", {
    set.seed(103)
    m <- matrix(runif(15), nrow = 5)
    cr <- constraints(m, ratio = "sequential")
    res <- recover_table(cr$p, colSums(m), rowSums(m), ratio = "sequential")
    expect_that(sum(abs(res$table-m)), is_less_than(1e-8))

})
