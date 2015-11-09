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
#'
#' @return a matrix
jac_cr <- function(x, A, B) {
    rbind(A,B/x)
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
constraints <- function(a) {
    A <- genA(nrow(a),ncol(a))
    B <- genB(nrow(a),ncol(a))
    cs <- colSums(a)
    rs <- rowSums(a)
    p <- B %*% log(as.vector(a))
    list(cs = cs, rs = rs[-nrow(a)], p = p, A = A[-nrow(A), ], B = B)
}

#' Recover table given column and row sums and initial guess
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
#' res <- recover_table(ma, colSums(m), rowSums(m))
#' res$table - m
recover_table <- function(initM, col_sums, row_sums, ...) {

    if (nrow(initM) != length(row_sums)) stop("The number of rows in initial matrix is not the same as in row totals")
    if (ncol(initM) != length(col_sums)) stop("The number of columns in initial matrix is not the same as in column totals")

    cr <- constraints(initM)

    o <- nleqslv(as.vector(initM), slv, jacobian = jac_cr,
            cs = col_sums, rs = row_sums[-length(row_sums)], p = cr$p,
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



