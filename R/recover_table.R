#' Generate a constraint matrix for column and row sums
#'
#' @param n integer, number of rows
#' @param m integer, number of columns
#'
#' @return constraint matrix with n+m rows and n*m columns
#' @import Matrix
genA <- function(n,m) {
    cind <- lapply((0:(m - 1))*n, function(shift) shift + 1:n)
    rind <- lapply(0:(n - 1), function(shift) shift + seq(1,(m - 1)*n + 1,by = n))
    rr1 <- rep(1:length(cind), times = sapply(cind,length))
    cc1 <- unlist(cind)
    rr2 <- rep(1:length(rind), times = sapply(rind,length)) + length(cind)
    cc2 <- unlist(rind)
    do.call("sparseMatrix",list(i = c(rr1,rr2),j = c(cc1,cc2), x = rep(1,length(rr1) + length(rr2))))
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
            resl[[i + (j - 1)*(n - 1)]] <- res
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
    rr <- rep(1:((n - 1)*(m - 1)),each = 4)
    cc <- rep(NA,length(rr))
    x <- rep(c(-1,1,1,-1),(n - 1)*(m - 1))
    for (j in 1:(m - 1)) {
        for (i in 1:(n - 1)) {
            cc[4*(i + (j - 1) * (n - 1) - 1) + 1] <- 1 #res[1] <- -1
            cc[4*(i + (j - 1) * (n - 1) - 1) + 2] <- i + 1# res[i + 1] <- 1
            cc[4*(i + (j - 1) * (n - 1) - 1) + 3] <- 1 + j*n # res[1 + j*n] <- 1
            cc[4*(i + (j - 1) * (n - 1) - 1) + 4] <- i + 1 + j*n
        }
    }
    do.call("sparseMatrix",list(i = rr,j = cc,x = x))
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

#' Recover table from a sequential log ratios and first column and row.
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


#' Recover table from a fixed log ratios and first column and row.
#'
#' @param p a matrix of log ratios
#' @param col1 a vector containing the first column
#' @param row1 a vector containing the first row
#'
#' @return a matrix
invp1 <- function(p, col1, row1) {
    #if (col1[1] != row1[1]) stop("The first element of first column should coincide with the first element of the first row")
    col1[1] <- row1[1]
    n <- length(col1)
    m <- length(row1)
    if (!identical(dim(p),as.integer(c(n - 1, m - 1)))) stop("The log ratio matrix dimensions are wrong")
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
#' @param p the matrix of fixed log ratios
#' @param exclude the matrix containing the indices of elements that should be excluded. Default is \code{NULL}, for no exclusions.
#'
#' @return a vector
cr <- function(x, p, exclude = NULL) {
    n <- nrow(p) + 1
    m <- ncol(p) + 1
    cc <- x[1:n]
    rr <- c(x[1], x[(n + 1):length(x)])
    M <- invp1(p, cc, rr)
    if (!is.null(exclude)) M[exclude] <- 0
    c(colSums(M),rowSums(M)[-n])
}

#' Jacobian of the system for recovering the table
#'
#' @param x the vectorized matrix
#' @param A the constraint matrix for sums
#' @param p the matrix of fixed log ratios
#' @param exclude the matrix containing the indices of elements that should be excluded. Default is \code{NULL}, for no exclusions.
#' @param ... additional parameters which are ignored.
#'
#' @return a matrix
jac_cr <- function(x, p,  A, exclude = NULL, ...) {
    n <- nrow(p) + 1
    m <- ncol(p) + 1

    res <- vector("list",n*m)
    cc <- x[1:n]
    rr <- c(x[1],x[(n + 1):(n + m - 1)])
    for (i in 1:n) {
        for (j in 1:m) {
            if (j == 1) {
                res[[i + (j - 1)*n]] <- list(j = i, x = 1)
            }
            if (i == 1 & j > 1) {
                res[[i + (j - 1)*n]] <- list(j = n + j - 1, x = 1)
            }
            if (i > 1 & j > 1) {
                dij <- exp(-p[i - 1, j - 1])
                res[[i + (j - 1)*n]] <- list(j = c(1, i, n + j - 1),
                                             x = dij/cc[1] *c(-cc[i] * rr[j]/ cc[1], rr[j], cc[i]))
            }
        }
    }
    ii <- rep(1:(n*m), times = sapply(res, function(x) length(x$j)))
    jj <- unlist(lapply(res,"[[", "j"))
    xx <- unlist(lapply(res,"[[", "x"))
    ff <- do.call("sparseMatrix",list(i = ii, j = jj, x = xx))
    #list(ff=ff,j=A %*% ff)

    if (!is.null(exclude)) {
        kk <- exclude[,1] + (exclude[,2] - 1)*n
        ff <- ff[-kk,]
        A <- A[, -kk]
    }
    as.matrix(A %*% ff)
}

#' The mismatch between the matrix and its constraints
#'
#' @param x the vectorized matrix
#' @param cs the vector of column sums
#' @param rs the vector of row sums
#' @param p the matrix of fixed log ratios
#' @param exclude the matrix containing the indices of elements that should be excluded. Default is \code{NULL}, for no exclusions.
#' @param ... additional parameters which are ignored
#'
#' @return the vector of differences
slv <- function(x, cs, rs, p, exclude = NULL, ...) {
    cr(x, p, exclude = exclude) - c(cs,rs)
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
    p <- as.vector(B %*% log(as.vector(a)))
    list(cs = cs, rs = rs[-nrow(a)], p = matrix(p, nrow = nrow(a) - 1), A = A[-nrow(A), ], ix = c(a[,1],a[1,-1]), ratio = ratio)

}

get_r <- function(cc, cs, rs, p, exclude = NULL) {
    d <- exp(-p)
    if (!is.null(exclude)) d[exclude] <- 0
    D <- rbind(rep(1,ncol(p)),d)
    as.vector(cs[-1]*cc[1]/(t(cc) %*% D))
}

get_c <- function(rr, cs, rs, p, exclude = NULL) {
    d <- exp(-p)
    if (!is.null(exclude)) d[exclude] <- 0
    rs <- c(rs, sum(cs) - sum(rs))
    c1 <- rs[1] - sum(rr)
    c2 <- c1*rs[-1]/(c1 + d %*% rr)
    c(c1,c2)
}

get_ij <- function(kk, n, m) {
    i <- kk %% n
    j <- kk %/% n + 1
    j[i == 0] <- j[i == 0] - 1
    i[i == 0] <- n
    cbind(i, j)
}

#' Recover table given column and row sums and given the log ratios
#'
#' @param p the log ratios
#' @param col_sums the vector of column sums
#' @param row_sums the vector of row sums
#' @param ratio the type of ratios used. The default is \code{"fixed"} for ratios with fixed first row and column.
#' It is also possible to use \code{"sequential"} ratios
#' @param exclude the matrix containing the indices of elements that should be excluded. Default is \code{NULL}, for no exclusions.
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
recover_table <- function(p, col_sums, row_sums, exclude = NULL, ratio = c("fixed", "sequential"), ...) {

    ratio <- match.arg(ratio)
    if (nrow(p) != length(row_sums) - 1 ) stop("The number of rows in log-likelihood ratio matrix should be one less than the number of row sums totals")
    if (ncol(p) != length(col_sums) - 1 ) stop("The number of columns in initial matrix log-likelihood ratio matrix should be one less than the number of column sums totals")
    if(ratio == "sequential" & !is.null(exclude)) stop("Table with fixed elements can only be recovered only fixed ratios")


    n <- nrow(p) + 1
    m <- ncol(p) + 1

    A <- genA(n, m)
    col1 <- c(sum(row_sums)/(n*m),row_sums[-1]/m)
    col1 <- col1*col_sums[1]/sum(col1)
    row1 <- get_r(col1, col_sums, row_sums[-n], p, exclude - 1)
    row1 <- row1*row_sums[1]/(col1[1] + sum(row1))
    col1 <- get_c(row1, col_sums, row_sums[-n], p, exclude - 1)
    col1 <- col1*col_sums[1]/sum(col1)

    ic <- c(sum(row_sums)/(n*m),col_sums[-1]/n)
    if (ratio == "sequential") p <- cumsum2(p)
    o <- nleqslv(c(col1, row1), slv, jac = jac_cr,
                 cs = col_sums, rs = row_sums[-length(row_sums)], p = p, A = A[-nrow(A), ], exclude = exclude, ...)
    if (o$termcd > 2) warning("The acceptable solution was not found. Numerical optimisation ended with the following message: ", o$message)

    ans <- invp1(p, o$x[1:n], c(o$x[1], o$x[(n + 1):(n + m - 1)]))
    if (!is.null(exclude)) ans[exclude] <- 0
    res <- list(table = ans, opt = o)
    class(res) <- "recover_table"
    res
}

#' Rescale table given column and row sums and initial table
#'
#' @param initM the initial matrix
#' @param col_sums the vector of column sums
#' @param row_sums the vector of row sums
#' @param exclude the matrix containing the indices of elements that should be excluded. Default is \code{NULL}, for no exclusions.
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
rescale_table <- function(initM, col_sums, row_sums, exclude = NULL, ...) {

    if (nrow(initM) != length(row_sums)) stop("The number of rows in initial matrix is not the same as in row totals")
    if (ncol(initM) != length(col_sums)) stop("The number of columns in initial matrix is not the same as in column totals")

    cr <- constraints(initM)
    nai <- which(is.na(cr$p) | is.infinite(cr$p))
    if (length(nai) > 0) {
        nae <- get_ij(nai, nrow(cr$p), ncol(cr$p)) + 1
        nae <- nae[order(nae[,1], nae[,2]), ]
        if (is.null(exclude)) stop("Please set exclusion for zero elements of the matrix")
        exclude <- exclude[order(exclude[,1], exclude[,2]), ]
        if (!identical(as.integer(exclude), as.integer(nae))) stop("Excluding restricting does not coincide with zero elements in the matrix")
    }
    o <- nleqslv(cr$ix, slv, jac = jac_cr,
                 cs = col_sums, rs = row_sums[-length(row_sums)], p = cr$p,
                 A = cr$A, exclude = exclude, ...)
    if (o$termcd > 2) warning("The acceptable solution was not found. Numerical optimisation ended with the following message: ", o$message)
    n <- nrow(initM)
    m <- ncol(initM)
    ans <- invp1(cr$p, o$x[1:n], c(o$x[1], o$x[(n + 1):(n + m - 1)]))
    if (!is.null(exclude)) ans[exclude] <- 0

    res <- list(table = ans, opt = o)
    class(res) <- "recover_table"
    res
}

#' @export
#' @method print recover_table
print.recover_table <- function(x, ...) {
    print(x$table)
}



