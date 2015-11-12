# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'



#' Generate a constraint matrix for column and row sums
#'
#' @param n integer, number of rows
#' @param m integer, number of columns
#'
#' @return constraint matrix with n+m rows and n*m columns
genA <- function(n,m) {
    cind <- lapply((0:(m - 1))*n, function(shift) shift + 1:n)
    rind <- lapply(0:(n - 1), function(shift) shift + seq(1,(m - 1)*n + 1,by = n))
    add_ones <- function(i) {
        res <- rep(0,n*m)
        res[i] <- 1
        res
    }
    do.call("rbind",lapply(c(cind,rind), add_ones))
}

#' Generate a constraint matrix for bivariate log-likelihood ratio
#'
#' @param n integer, number of rows
#' @param m integer, number of columns
#'
#' @return constraint matrix with (n-1)*(m-1) rows and n*m columns
genB <- function(n,m) {
    resl <- vector("list",(n - 1)*(m - 1))
    for (i in 1:(n - 1)) {
        for (j in 1:(m - 1)) {
            res <- rep(0,n*m)
            res[i + (j - 1)*n] <- -1
            res[i + 1 + (j - 1)*n] <- 1
            res[i + j*n] <- 1
            res[i + 1 + j*n] <- -1
            resl[[i + (j - 1)*n]] <- res
        }
    }
    do.call("rbind",resl)
}

#' Generate a constraint matrix with fixed first row and column
#' for bivariate log-likelihood ratio
#'
#' @param n integer, number of rows
#' @param m integer, number of columns
#'
#' @return constraint matrix with (n-1)*(m-1) rows and n*m columns
genB1 <- function(n,m) {
    resl <- vector("list",(n - 1)*(m - 1))
    for (i in 1:(n - 1)) {
        for (j in 1:(m - 1)) {
            res <- rep(0,n*m)
            res[1] <- -1
            res[i + 1] <- 1
            res[1 + j*n] <- 1
            res[i + 1 + j*n] <- -1
            resl[[i + (j - 1)*n]] <- res
        }
    }
    do.call("rbind",resl)
}

#' Cumulative sum for a matrix
#'
#' If p is initial matrix, the result is the matrix with the following elements:
#' \deqn{r_{ij} = \sum_{k=1}^i\sum_{l=1}^j p_{ij}}
#'
#' @param p a matrix
#'
#' @return a matrix of cumulative sums
cumsum2 <- function(p) {
    n <- nrow(p)
    m <- ncol(p)
    res <- matrix(0,nrow = n, ncol = m)
    for (j in 1:m) {
        for (i in 1:n) {
            res[i, j] <- sum(p[1:i,1:j])
        }
    }
    res
}

#' Recover table from a sequential log-likelihood ratios and first column and row.
#'
#' @param p a matrix of log-likelihood ratios
#' @param col1 a vector containing the first column
#' @param row1 a vector containing the first row
#'
#' @return a matrix
invp <- function(p, col1, row1) {
    if (col1[1] != row1[1]) stop("The first element of first column should coincide with the first element of the first row")
    n <- length(col1)
    m <- length(row1)
    if (!identical(dim(p),as.integer(c(n - 1, m - 1)))) stop("The log likelihood ratio matrix dimensions are wrong")
    res <- matrix(0, ncol = m, nrow = n)
    res[1,] <- row1
    res[,1] <- col1
    cc <- outer(log(col1), log(row1), "+")[2:n,2:m]
    res[2:n,2:m] <- exp(cc - log(row1[1]) - cumsum2(p))
    res
}

#' Recover table from a fixed log-likelihood ratios and first column and row.
#'
#' @param p a matrix of log-likelihood ratios
#' @param col1 a vector containing the first column
#' @param row1 a vector containing the first row
#'
#' @return a matrix
invp1 <- function(p, col1, row1) {
    if (col1[1] != row1[1]) stop("The first element of first column should coincide with the first element of the first row")
    n <- length(col1)
    m <- length(row1)
    if (!identical(dim(p),as.integer(c(n - 1, m - 1)))) stop("The log likelihood ratio matrix dimensions are wrong")
    res <- matrix(0, ncol = m, nrow = n)
    res[1,] <- row1
    res[,1] <- col1
    cc <- outer(log(col1), log(row1), "+")[2:n,2:m]
    res[2:n,2:m] <- exp(cc - log(row1[1]) - p)
    res
}


#' Generate the constraints given the vectorized matrix and the constraint matrices
#'
#' @param x the vectorized matrix
#' @param A the constraint matrix for sums
#' @param B the constraint matrix for log-likelihood ratio
#'
#' @return a vector
cr <- function(x,A,B) {
    c(A %*% x, B %*% log(x))
}

#' Jacobian of the system for recovering the table
#'
#' @param x the vectorized matrix
#' @param A the constraint matrix for sums
#' @param B the constraint matrix for log-likelihood ratio
#' @param ... additional parameters which are ignored.
#'
#' @return a matrix
jac_cr <- function(x, A, B, ...) {
    rbind(A,t(t(B)/x))
}

#' The mismatch between the matrix and its constraints
#'
#' @param x the vectorized matrix
#' @param cs the vector of column sums
#' @param rs the vector of row sums
#' @param p the vector of log-likelihood ratios
#' @param A the constraint matrix for sums
#' @param B the constraint matrix for log-likelihood ratios
#'
#' @return the vector of differences
slv <- function(x, cs, rs, p, A, B) {
    cr(x,A,B) - c(cs,rs,p)
}

#' Calculate the matrix constraints
#'
#' @param a the matrix
#' @param ratio a string specifying the type of ratios used. The default is \code{"fixed"} for ratios with fixed first row and column.
#' It is also possible to use \code{"sequential"} ratios
#'
#' @return a list with the following elements
#' \item{cs}{the vector of column sums}
#' \item{rs}{the vector of row sums without the last element}
#' \item{p}{the vector of likelihood ratios}
#' \item{A}{the constraint matrix for sums without the last row}
#' \item{B}{the constraint matrix for log-likelihood ratios}
#' @export
#' @examples
#' m <- matrix(1:20, nrow = 5)
#' constraints(m)
constraints <- function(a, ratio = c("fixed", "sequential")) {
    ratio <- match.arg(ratio)
    A <- genA(nrow(a),ncol(a))
    B <- switch(ratio,
                fixed = genB1(nrow(a),ncol(a)),
                sequential = genB(nrow(a),ncol(a)))

    cs <- colSums(a)
    rs <- rowSums(a)
    p <- B %*% log(as.vector(a))
    list(cs = cs, rs = rs[-nrow(a)], p = matrix(p, nrow = nrow(a) - 1), A = A[-nrow(A), ], B = B, ratio = ratio)
}

#' Recover table given column and row sums and given the log ratios
#'
#' @param p the log ratios
#' @param col_sums the vector of column sums
#' @param row_sums the vector of row sums
#' @param ratio the type of ratios used. The default is \code{"fixed"} for ratios with fixed first row and column.
#' It is also possible to use \code{"sequential"} ratios
#' @param ... additional arguments to nleqlsv
#'
#' @return a list with the following elements
#' \item{table}{a recovered table}
#' \item{opt}{a result of \code{nleqslv}}
#' @export
#' @import nleqslv
#' @examples
#' set.seed(10)
#' m <- matrix(runif(9),3)
#' cr <- constraints(m)
#' res <- recover_table(cr$p, colSums(m), rowSums(m))
#' res$table - m
recover_table <- function(p, col_sums, row_sums, ratio = c("fixed", "sequential"), ...) {

    ratio <- match.arg(ratio)
    if (nrow(p) != length(row_sums) - 1 ) stop("The number of rows in log-likelihood ratio matrix should be one less than the number of row sums totals")
    if (ncol(p) != length(col_sums) - 1 ) stop("The number of columns in initial matrix log-likelihood ratio matrix should be one less than the number of column sums totals")

    n <- nrow(p) + 1
    m <- ncol(p) + 1
    A <- genA(n, m)
    B <- switch(ratio,
                fixed = genB1(n,m),
                sequential = genB(n,m))

    ir <- c(sum(row_sums)/(n*m),row_sums[-1]/m)
    ic <- c(sum(row_sums)/(n*m),col_sums[-1]/n)
    initM <- switch(ratio,
                    fixed = invp1(p, ir, ic),
                    sequential = invp1(p, ir, ic)
    )
    o <- nleqslv(as.vector(initM), slv, jac = jac_cr,
            cs = col_sums, rs = row_sums[-length(row_sums)], p,
            A = A[-nrow(A),], B = B, ...)
    if (o$termcd > 2) warning("The acceptable solution was not found. Numerical optimisation ended with the following message: ", o$message)
    ans <- matrix(o$x, nrow = n)
    res <- list(table = ans, opt = o)
    class(res) <- "recover_table"
    res
}


#' Rescale table given column and row sums and initial table
#'
#' @param initM the initial matrix
#' @param col_sums the vector of column sums
#' @param row_sums the vector of row sums
#' @param ... additional arguments to nleqlsv
#'
#' @return a list with the following elements
#' \item{table}{a recovered table}
#' \item{opt}{a result of \code{nleqslv}}
#' @export
#' @import nleqslv
#' @examples
#' set.seed(10)
#' m <- matrix(runif(9),3)
#' ma <- m + matrix(runif(9),3)/20
#' res <- rescale_table(ma, colSums(m), rowSums(m))
#' res$table - m
rescale_table <- function(initM, col_sums, row_sums, ...) {

    if (nrow(initM) != length(row_sums)) stop("The number of rows in initial matrix is not the same as in row totals")
    if (ncol(initM) != length(col_sums)) stop("The number of columns in initial matrix is not the same as in column totals")

    cr <- constraints(initM)

    o <- nleqslv(as.vector(initM), slv, jac = jac_cr,
                 cs = col_sums, rs = row_sums[-length(row_sums)], p = as.vector(cr$p),
                 A = cr$A, B = cr$B, ...)
    if (o$termcd > 2) warning("The acceptable solution was not found. Numerical optimisation ended with the following message: ", o$message)
    ans <- matrix(o$x, nrow = nrow(initM))
    res <- list(table = ans, opt = o)
    class(res) <- "recover_table"
    res
}

#' @export
#' @method print recover_table
print.recover_table <- function(x, ...) {
    print(x$table)
}



