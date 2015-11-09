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

cr <- function(x,A,B) {
    c(A %*% x, B %*% log(x))
}

jac_cr <- function(x, A, B) {
    rbind(A,B/x)
}

slv <- function(x, cs, rs, p, A, B) {
    cr(x,A,B) - c(cs,rs,p)
}

constraints <- function(a) {
    A <- genA(nrow(a),ncol(a))
    B <- genB(nrow(a),ncol(a))
    cs <- colSums(a)
    rs <- rowSums(a)
    p <- B %*% log(as.vector(a))
    list(cs = cs, rs = rs[-nrow(a)], p = p, A = A[-nrow(A), ], B = B)
}

recover_table <- function(initM, col_sums, row_sums, ...) {

    if (nrow(initM) != length(row_sums)) stop("The number of rows in initial matrix is not the same as in row totals")
    if (ncol(initM) != length(col_sums)) stop("The number of columns in initial matrix is not the same as in column totals")

    cr <- constraints(initM)

    o <- nleqslv(as.vector(initM), slv, jacobian = jac_cr,
            cs = col_sums, rs = row_sums[-length(row_sums)], p = cr$p,
            A = cr_m$A, B = cr_m$B, ...)
    if (o$termcmd > 2) warning("The acceptable solution was not found. Numerical optimisation ended with the following message: ", o$message)
    ans <- matrix(o$x, nrow = nrow(initM))
    res <- list(table = ans, opt = o)
    class(res) <- "recover_table"
}

print.recover_table <- function(x, ...) {
    print(x$table)
}



