context("Testing exclusions")

sample_ij <- function(ni, n, m) {
    fc <- 1:n
    fr <- (1:(m - 1))*n + 1
    kk <- sample((1:(n*m))[-c(fc, fr)], ni, replace = FALSE )
    i <- kk %% n
    j <- kk %/% n + 1
    j[i == 0] <- j[i == 0] - 1
    i[i == 0] <- n
    if (1 %in% i | 1 %in% j ) stop("You cannot have zero elements in first column and row")
    cbind(i, j)
}


test_that("The jacobian works", {
    set.seed(1001)
    m <- matrix(runif(63),nrow = 9)

    exc <- sample_ij(5, 9, 7)
    m[exc] <- 0
    cr <- constraints(m)
    jc_t <- jac_cr(cr$ix, p = cr$p, A = cr$A, exclude = exc)
    jc_n <- jacobian(slv, cr$ix, cs = cr$cs, rs = cr$rs, p = cr$p, exclude = exc)
    expect_that( sum(abs(jc_t - jc_n)), is_less_than(1e-6))
})

test_that("Recovering the matrix with exclusions given the constraints works", {
    set.seed(1002)
    m <- matrix(runif(15), nrow = 5)
    exc <- sample_ij(6, 5, 3)
    m[exc] <- 0

    cr <- constraints(m)
    o <- nleqslv(cr$ix + runif(length(cr$ix))/14, slv, jacobian = jac_cr,
                 cs = cr$cs, rs = cr$rs, p = cr$p, A = cr$A, exclude = exc)

    res <- invp1(cr$p, o$x[1:5], c(o$x[1], o$x[6:7]))
    res[exc] <- 0
    expect_that(sum(abs(res - m)), is_less_than(1e-6))
})


test_that("Trivial rescaling with exclusions works", {
    set.seed(1003)
    m <- matrix(runif(24), nrow = 4)
    exc <- sample_ij(8, 4, 6)
    m[exc] <- 0
    res <- rescale_table(m, colSums(m), rowSums(m), exclude = exc)
    expect_that(sum(abs(res$table - m)), is_less_than(1e-8))
})

test_that("Trivial recovering works with exclusions works", {
    set.seed(1004)
    mm <- matrix(runif(36), nrow = 3)
    exc <- sample_ij(10, 3, 12)
    mm[exc] <- 0
    cr <-  constraints(mm, ratio = "fixed")
    res <- recover_table(cr$p, colSums(mm), rowSums(mm), ratio = "fixed", exclude = exc )
    expect_that(sum(abs(res$table - mm)), is_less_than(1e-6))
})


