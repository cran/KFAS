estmiss <-
function (out) 
{
    Ft <- out$Ft
    Fstar <- out$Fstar
    Finf <- out$Finf
    tv <- NULL
    tv[1] <- !(is.na(dim(as.array(Zt))[3]) || dim(as.array(Zt))[3] == 
        1)
    tv[2] <- !(is.na(dim(as.array(Ht))[3]) || dim(as.array(Ht))[3] == 
        1)
    Zt <- array(out$Zt, c(out$p, out$m, (out$n - 1) * tv[1] + 
        1))
    Ht <- array(out$Ht, c(out$p, out$p, (out$n - 1) * tv[2] + 
        1))
    if (out$d > 0) {
        for (i in 1:out$d) {
            if (out$ydimt[i] != out$p) {
                Fstar[, , i] <- as.matrix(Zt[, , (i - 1) * tv[1] + 
                  1]) %*% out$Pstar[, , i] %*% t(as.matrix(Zt[, 
                  , (i - 1) * tv[1] + 1])) + Ht[, , (i - 1) * 
                  tv[2] + 1]
                Finf[, , i] <- as.matrix(Zt[, , (i - 1) * tv[1] + 
                  1]) %*% out$Pinf[, , i] %*% t(as.matrix(Zt[, 
                  , (i - 1) * tv[1] + 1]))
            }
        }
    }
    if (out$d < out$n) {
        for (i in (out$d + 1):out$n) {
            if (out$ydimt[i] != out$p) 
                Ft[, , i] <- as.matrix(Zt[, , (i - 1) * tv[1] + 
                  1]) %*% out$Pt[, , i] %*% t(as.matrix(Zt[, 
                  , (i - 1) * tv[1] + 1])) + Ht[, , (i - 1) * 
                  tv[2] + 1]
        }
    }
    if (!is.null(out$ahat)) {
        yt <- out$yt
        for (i in 1:out$n) {
            if (out$ydimt[i] != out$p) 
                yt[, i] <- as.matrix(Zt[, , (i - 1) * tv[1] + 
                  1]) %*% as.matrix(out$ahat[, i])
        }
        est.out <- list(Fstar = Fstar, Finf = Finf, Ft = Ft, 
            yt = yt)
    }
    else {
        est.out <- list(Fstar = Fstar, Finf = Finf, Ft = Ft)
    }
    return(est.out)
}
forecast <-
function (out, fc = 1, Zt.fc = NULL, Tt.fc = NULL, Rt.fc = NULL, 
    Ht.fc = NULL, Qt.fc = NULL) 
{
    yt.fc <- array(0, dim = c(out$p, fc))
    Ft.fc <- array(0, dim = c(out$p, out$p, fc))
    at.fc <- array(0, dim = c(out$m, fc + 1))
    Pt.fc <- array(0, dim = c(out$m, out$m, fc + 1))
    at.fc[, 1] <- out$at[, out$n + 1]
    Pt.fc[, , 1] <- out$Pt[, , out$n + 1]
    tv <- rep(1, 5)
    error <- ""
    if (is.null(Zt.fc)) {
        if (!(is.na(dim(as.array(out$Zt))[3]) || dim(as.array(out$Zt))[3] == 
            1)) {
            error <- paste(error, "Zt is not time-invariant and Zt.fc is NULL!")
        }
        else {
            Zt.fc <- array(out$Zt, c(out$p, out$m, 1))
            tv[1] <- 0
        }
    }
    if (is.null(Tt.fc)) {
        if (!(is.na(dim(as.array(out$Tt))[3]) || dim(as.array(out$Tt))[3] == 
            1)) {
            error <- paste(error, "Tt is not time-invariant and Tt.fc is NULL!")
        }
        else {
            Tt.fc <- array(out$Tt, c(out$m, out$m, 1))
            tv[2] <- 0
        }
    }
    if (is.null(Rt.fc)) {
        if (!(is.na(dim(as.array(out$Rt))[3]) || dim(as.array(out$Rt))[3] == 
            1)) {
            error <- paste(error, "Rt is not time-invariant and Rt.fc is NULL!")
        }
        else {
            Rt.fc <- array(out$Rt, c(out$m, out$r, 1))
            tv[3] <- 0
        }
    }
    if (is.null(Ht.fc)) {
        if (!(is.na(dim(as.array(out$Ht))[3]) || dim(as.array(out$Ht))[3] == 
            1)) {
            error <- paste(error, "Ht is not time-invariant, and Ht.fc is NULL!")
        }
        else {
            Ht.fc <- array(out$Ht, c(out$p, out$p, 1))
            tv[4] <- 0
        }
    }
    if (is.null(Qt.fc)) {
        if (!(is.na(dim(as.array(out$Qt))[3]) || dim(as.array(out$Qt))[3] == 
            1)) {
            error <- paste(error, "Qt is not time-invariant, and Qt.fc is NULL!")
        }
        else {
            Qt.fc <- array(out$Qt, c(out$r, out$r, 1))
            tv[5] <- 0
        }
    }
    if (error != "") 
        stop(error)
    for (i in 1:fc) {
        yt.fc[, i] <- matrix(Zt.fc[, , ((i - 1) * tv[1] + 1)], 
            out$p, out$m) %*% matrix(at.fc[, i], out$m, 1)
        Ft.fc[, , i] <- matrix(Zt.fc[, , ((i - 1) * tv[1] + 1)], 
            out$p, out$m) %*% Pt.fc[, , i] %*% t(matrix(Zt.fc[, 
            , ((i - 1) * tv[1] + 1)], out$p, out$m)) + Ht.fc[, 
            , ((i - 1) * tv[4] + 1)]
        at.fc[, i + 1] <- Tt.fc[, , ((i - 1) * tv[2] + 1)] %*% 
            matrix(at.fc[, i], out$m, 1)
        Pt.fc[, , i + 1] <- Tt.fc[, , ((i - 1) * tv[2] + 1)] %*% 
            Pt.fc[, , i] %*% t(Tt.fc[, , ((i - 1) * tv[2] + 1)]) + 
            matrix(Rt.fc[, , ((i - 1) * tv[3] + 1)], out$m, out$r) %*% 
                Qt.fc[, , ((i - 1) * tv[5] + 1)] %*% t(matrix(Rt.fc[, 
                , ((i - 1) * tv[3] + 1)], out$m, out$r))
    }
    return(list(yt.fc = yt.fc, Ft.fc = Ft.fc, at.fc = at.fc, 
        Pt.fc = Pt.fc))
}
kf <-
function (yt, Zt, Tt, Rt, Ht, Qt, a1, P1, P1inf = 0, optcal = c(TRUE, 
    TRUE, TRUE, TRUE), tol = 1e-07) 
{
    if (!is.array(yt)) {
        if (!is.matrix(yt)) 
            yt <- array(yt, dim = c(1, length(yt)))
        else yt <- array(yt, dim = dim(yt))
    }
    storage.mode(yt) <- "double"
    d <- 0
    p <- dim(yt)[1]
    n <- dim(yt)[2]
    m <- length(a1)
    if (is.vector(Qt)) 
        r <- 1
    else r <- dim(as.array(Qt))[2]
    tv <- array(0, dim = 3)
    tv[1] <- !(is.na(dim(as.array(Tt))[3]) || dim(as.array(Tt))[3] == 
        1)
    tv[2] <- !(is.na(dim(as.array(Rt))[3]) || dim(as.array(Rt))[3] == 
        1)
    tv[3] <- !(is.na(dim(as.array(Qt))[3]) || dim(as.array(Qt))[3] == 
        1)
    ymiss <- is.na(yt)
    ydimt <- array(0, dim = n)
    H <- array(Ht, c(p, p, n))
    Z <- array(Zt, dim = c(p, m, n))
    y <- yt
    for (i in 1:n) {
        ydimt[i] <- sum(!ymiss[1:p, i])
        if (ydimt[i] != p && ydimt[i] != 0) {
            y[1:ydimt[i], i] <- yt[!ymiss[, i], i]
            H[1:ydimt[i], 1:ydimt[i], i] <- H[!ymiss[, i], !ymiss[, 
                i], i]
            Z[1:ydimt[i], , i] <- Z[!ymiss[, i], , i]
        }
    }
    at <- array(0, dim = c(m, n + 1))
    Pt <- array(0, dim = c(m, m, n + 1))
    vt <- array(0, dim = c(p, n))
    vtuni <- array(0, dim = c(p, n))
    Ft <- array(0, dim = c(p, p, n))
    Ftuni <- array(0, dim = c(p, n))
    Kt <- array(0, dim = c(m, p, n))
    Ktuni <- array(0, dim = c(m, p, n))
    Lt <- array(0, dim = c(m, m, n))
    Pinf <- array(0, dim = c(m, m, n + 1))
    Pstar <- array(0, dim = c(m, m, n + 1))
    Kinf <- array(0, dim = c(m, p, n))
    Kstar <- array(0, dim = c(m, p, n))
    Kinfuni <- array(0, dim = c(m, p, n))
    Kstaruni <- array(0, dim = c(m, p, n))
    Finfuni <- array(0, dim = c(p, n))
    Fstaruni <- array(0, dim = c(p, n))
    Finf <- array(0, dim = c(p, p, n))
    Fstar <- array(0, dim = c(p, p, n))
    Linf <- array(0, dim = c(m, m, n))
    Lstar <- array(0, dim = c(m, m, n))
    Pinf[, , 1] <- P1inf
    storage.mode(Z) <- "double"
    storage.mode(H) <- "double"
    storage.mode(y) <- "double"
    storage.mode(d) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(r) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(tv) <- "integer"
    storage.mode(Zt) <- "double"
    storage.mode(Tt) <- "double"
    storage.mode(Rt) <- "double"
    storage.mode(Ht) <- "double"
    storage.mode(Qt) <- "double"
    storage.mode(a1) <- "double"
    storage.mode(P1) <- "double"
    lik <- 0
    storage.mode(lik) <- "double"
    info <- 0
    storage.mode(info) <- "integer"
    storage.mode(optcal) <- "integer"
    storage.mode(ydimt) <- "integer"
    storage.mode(vt) <- "double"
    storage.mode(vtuni) <- "double"
    storage.mode(Ft) <- "double"
    storage.mode(Ftuni) <- "double"
    storage.mode(Finfuni) <- "double"
    storage.mode(Fstaruni) <- "double"
    storage.mode(Finf) <- "double"
    storage.mode(Fstar) <- "double"
    storage.mode(Kt) <- "double"
    storage.mode(Ktuni) <- "double"
    storage.mode(Kinfuni) <- "double"
    storage.mode(Kinf) <- "double"
    storage.mode(Kstaruni) <- "double"
    storage.mode(Kstar) <- "double"
    storage.mode(Pinf) <- "double"
    storage.mode(Pstar) <- "double"
    kfout <- NULL
    j <- 0
    storage.mode(j) <- "integer"
    kfout <- .Fortran("kf", PACKAGE = "KFAS", NAOK = TRUE, yt = y, 
        ydimt = ydimt, tv = tv, Zt = Z, Tt = array(Tt, c(m, m, 
            (n - 1) * tv[1] + 1)), Rt = array(Rt, c(m, r, (n - 
            1) * tv[2] + 1)), Ht = H, Qt = array(Qt, c(r, r, 
            (n - 1) * tv[3] + 1)), a1 = array(a1, c(m)), P1 = array(P1, 
            c(m, m)), at = at, Pt = Pt, vtuni = vtuni, Ftuni = Ftuni, 
        Ktuni = Ktuni, Pinf = Pinf, Pstar = Pstar, Finfuni = Finfuni, 
        Fstaruni = Fstaruni, Kinfuni = Kinfuni, Kstaruni = Kstaruni, 
        d = d, j = j, p = p, m = m, r = r, n = n, lik = lik, 
        optcal = optcal, info = info, vt = vt, Ft = Ft, Kt = Kt, 
        Lt = Lt, Finf = Finf, Fstar = Fstar, Kinf = Kinf, Kstar = Kstar, 
        Linf = Linf, Lstar = Lstar, tol = tol)
    kfout$Pinf <- array(kfout$Pinf[, , 1:(kfout$d + 1)], c(m, 
        m, (kfout$d + 1) * (kfout$d > 0)))
    kfout$Pstar <- array(kfout$Pstar[, , 1:(kfout$d + 1)], c(m, 
        m, (kfout$d + 1) * (kfout$d > 0)))
    kfout$Finfuni <- array(kfout$Finfuni[, 1:kfout$d], c(p, kfout$d))
    kfout$Fstaruni <- array(kfout$Fstaruni[, 1:kfout$d], c(p, 
        kfout$d))
    kfout$Kinfuni <- array(kfout$Kinfuni[, , 1:kfout$d], c(m, 
        p, kfout$d))
    kfout$Kstaruni <- array(kfout$Kstaruni[, , 1:kfout$d], c(m, 
        p, kfout$d))
    kfout$Finf <- array(kfout$Finf[, , 1:kfout$d], c(p, p, kfout$d))
    kfout$Fstar <- array(kfout$Fstar[, , 1:kfout$d], c(p, p, 
        kfout$d))
    kfout$Kinf <- array(kfout$Kinf[, , 1:kfout$d], c(m, p, kfout$d))
    kfout$Kstar <- array(kfout$Kstar[, , 1:kfout$d], c(m, p, 
        kfout$d))
    kfout$Linf <- array(kfout$Linf[, , 1:kfout$d], c(m, m, kfout$d))
    kfout$Lstar <- array(kfout$Lstar[, , 1:kfout$d], c(m, m, 
        kfout$d))
    kfout$yt <- yt
    kfout$Tt <- Tt
    kfout$Rt <- Rt
    kfout$Qt <- Qt
    kfout$Zt <- Zt
    kfout$Ht <- Ht
    if (kfout$d > 0) {
        kfout$Ft[, , 1:kfout$d] <- NA
        kfout$Ftuni[1:kfout$j, 1:kfout$d] <- NA
        kfout$Pt[, , 1:kfout$d] <- NA
        for (i in 1:kfout$d) {
            if (ydimt[i] != p) {
                kfout$Finfuni[(ydimt[i] + 1):p, i] <- NA
                kfout$Fstaruni[(ydimt[i] + 1):p, i] <- NA
                kfout$vtuni[(ydimt[i] + 1):p, i] <- NA
                kfout$Ktuni[, (ydimt[i] + 1):p, i] <- NA
                kfout$Kinfuni[, (ydimt[i] + 1):p, i] <- NA
                kfout$Kstaruni[, (ydimt[i] + 1):p, i] <- NA
                if (optcal[1]) 
                  kfout$vt[(ydimt[i] + 1):p, i] <- NA
                if (optcal[2]) {
                  kfout$Finf[(ydimt[i] + 1):p, (ydimt[i] + 1):p, 
                    i] <- NA
                  kfout$Fstar[(ydimt[i] + 1):p, (ydimt[i] + 1):p, 
                    i] <- NA
                }
                if (optcal[3] && optcal[2]) {
                  kfout$Kinf[, (ydimt[i] + 1):p, i] <- NA
                  kfout$Kstar[, (ydimt[i] + 1):p, i] <- NA
                }
            }
        }
    }
    if (kfout$d < n) {
        for (i in (kfout$d + 1):n) {
            if (ydimt[i] != p) {
                kfout$Ftuni[(ydimt[i] + 1):p, i] <- NA
                kfout$vtuni[(ydimt[i] + 1):p, i] <- NA
                kfout$Ktuni[, (ydimt[i] + 1):p, i] <- NA
                if (optcal[1]) 
                  kfout$vt[(ydimt[i] + 1):p, i] <- NA
                if (optcal[2]) 
                  kfout$Ft[(ydimt[i] + 1):p, (ydimt[i] + 1):p, 
                    i] <- NA
                if (optcal[3] && optcal[2]) 
                  kfout$Kt[, (ydimt[i] + 1):p, i] <- NA
            }
        }
    }
    if (optcal[1] == 0) 
        kfout$vt <- NULL
    if (optcal[2] == 0) {
        kfout$Ft <- NULL
        kfout$Finf <- NULL
        kfout$Fstar <- NULL
    }
    if (optcal[3] == 0) {
        kfout$Kt <- NULL
        kfout$Kinf <- NULL
        kfout$Kstar <- NULL
    }
    if (optcal[4] == 0) {
        kfout$Lt <- NULL
        kfout$Linf <- NULL
        kfout$Lstar <- NULL
    }
    if (kfout$info != 0) {
        if (kfout$info == 1) {
            kfout$lik <- -Inf
            warning("Could not diagonalize Ht")
        }
        else if (kfout$info == 2) {
            warning("Could not compute multivariate Kstar because multivariate Fstar is singular")
        }
        else if (kfout$info == 3) {
            warning("Could not compute multivariate Kinf because multivariate Finf is singular")
        }
        else warning("Could not compute multivariate Kt because multivariate Ft is singular")
    }
    return(kfout)
}
ks <-
function (out) 
{
    ahat <- array(0, dim = c(out$m, out$n))
    Vt <- array(0, dim = c(out$m, out$m, out$n))
    Nt <- array(0, dim = c(out$m, out$m, out$n + 1))
    Nt0 <- array(0, dim = c(out$m, out$m, out$d + 1))
    Nt1 <- array(0, dim = c(out$m, out$m, out$d + 1))
    Nt2 <- array(0, dim = c(out$m, out$m, out$d + 1))
    epshat <- array(0, dim = c(out$p, out$n))
    epshatvar <- array(0, dim = c(out$p, out$n))
    etahat <- array(0, dim = c(out$r, out$n))
    etahatvar <- array(0, dim = c(out$r, out$r, out$n))
    rt <- array(0, dim = c(out$m, out$n + 1))
    rt0 <- array(0, dim = c(out$m, out$d + 1))
    rt1 <- array(0, dim = c(out$m, out$d + 1))
    storage.mode(Nt0) <- "double"
    storage.mode(Nt1) <- "double"
    storage.mode(Nt2) <- "double"
    storage.mode(rt0) <- "double"
    storage.mode(rt1) <- "double"
    storage.mode(ahat) <- "double"
    storage.mode(Vt) <- "double"
    storage.mode(Nt) <- "double"
    storage.mode(rt) <- "double"
    ymiss <- is.na(out$yt)
    H <- array(out$Ht, c(out$p, out$p, out$n))
    Z <- array(out$Zt, dim = c(out$p, out$m, out$n))
    for (i in 1:out$n) {
        if (out$ydimt[i] != out$p && out$ydimt[i] != 0) {
            H[1:out$ydimt[i], 1:out$ydimt[i], i] <- H[!ymiss[, 
                i], !ymiss[, i], i]
            Z[1:out$ydimt[i], , i] <- Z[!ymiss[, i], , i]
        }
    }
    storage.mode(H) <- "double"
    storage.mode(Z) <- "double"
    ks.out <- .Fortran("ks", PACKAGE = "KFAS", NAOK = TRUE, out$ydimt, 
        out$tv, Z = Z, array(out$Tt, c(out$m, out$m, (out$n - 
            1) * out$tv[1] + 1)), H = H, array(out$Rt, c(out$m, 
            out$r, (out$n - 1) * out$tv[2] + 1)), array(out$Qt, 
            c(out$r, out$r, (out$n - 1) * out$tv[3] + 1)), out$at, 
        out$Pt, out$vtuni, out$Ftuni, out$Ktuni, ahat = ahat, 
        Vt = Vt, rt = rt, rt0 = rt0, rt1 = rt1, Nt = Nt, Nt0 = Nt0, 
        Nt1 = Nt1, Nt2 = Nt2, epshat = epshat, epshatvar = epshatvar, 
        etahat = etahat, etahatvar = etahatvar, Pinf = out$Pinf, 
        Pstar = out$Pstar, Kinfuni = out$Kinfuni, Kstaruni = out$Kstaruni, 
        Finfuni = out$Finfuni, Fstaruni = out$Fstaruni, d = out$d, 
        j = out$j, p = out$p, m = out$m, r = out$r, n = out$n, 
        tol = out$tol)
    ks.out <- c(out, list(ahat = ks.out$ahat, Vt = ks.out$Vt, 
        rt = ks.out$rt, rt0 = ks.out$rt0, rt1 = ks.out$rt1, Nt = ks.out$Nt, 
        Nt0 = ks.out$Nt0, Nt1 = ks.out$Nt1, Nt2 = ks.out$Nt2, 
        epshat = ks.out$epshat, epshatvar = ks.out$epshatvar, 
        etahat = ks.out$etahat, etahatvar = ks.out$etahatvar))
    ks.out
}
simsmoother <-
function (yt, Zt, Tt, Rt, Ht, Qt, a1, P1, P1inf = 0, nsim = 1, 
    tol = 1e-07) 
{
    if (!is.array(yt)) {
        if (!is.matrix(yt)) 
            cat("yt must be array or matrix")
        else yt <- array(yt, dim = dim(yt))
    }
    storage.mode(yt) <- "double"
    p <- dim(yt)[1]
    n <- dim(yt)[2]
    m <- length(a1)
    if (is.vector(Qt)) 
        r <- 1
    else r <- dim(as.array(Qt))[2]
    tv <- array(0, dim = 4)
    tv[1] <- !(is.na(dim(as.array(Tt))[3]) || dim(as.array(Tt))[3] == 
        1)
    tv[2] <- !(is.na(dim(as.array(Rt))[3]) || dim(as.array(Rt))[3] == 
        1)
    tv[3] <- !(is.na(dim(as.array(Qt))[3]) || dim(as.array(Qt))[3] == 
        1)
    ymiss <- is.na(yt)
    ydimt <- array(0, dim = n)
    H <- array(Ht, c(p, p, n))
    Z <- array(Zt, dim = c(p, m, n))
    y <- yt
    for (i in 1:n) {
        ydimt[i] <- sum(!ymiss[1:p, i])
        if (ydimt[i] != p && ydimt[i] != 0) {
            y[1:ydimt[i], i] <- yt[!ymiss[, i], i]
            H[1:ydimt[i], 1:ydimt[i], i] <- H[!ymiss[, i], !ymiss[, 
                i], i]
            Z[1:ydimt[i], , i] <- Z[!ymiss[, i], , i]
        }
    }
    Tt <- array(Tt, c(m, m, (n - 1) * tv[1] + 1))
    Rt <- array(Rt, c(m, r, (n - 1) * tv[2] + 1))
    Qt <- array(Qt, c(r, r, (n - 1) * tv[3] + 1))
    P1 <- array(P1, c(m, m))
    P1inf <- array(P1inf, c(m, m))
    a1 <- array(a1, c(m, 1))
    alphasim <- array(0, c(m, n, nsim))
    info <- 0
    epsplus <- array(0, c(p, n, nsim))
    etaplus <- array(0, c(r, n, nsim))
    aplus1 <- array(0, dim = c(m, nsim))
    storage.mode(tv) <- "integer"
    storage.mode(ydimt) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(r) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(info) <- "integer"
    storage.mode(nsim) <- "integer"
    for (t in 1:n) {
        if (ydimt[t] > 0) {
            epsplus[1:ydimt[t], t, ] <- rnorm(ydimt[t] * nsim, 
                mean = 0, sd = 1)
        }
    }
    etaplus[, , ] <- rnorm(r * n * nsim, mean = 0, sd = 1)
    nde<-which(diag(P1inf)==0)
    nnd<-length(nde)
    if (nnd > 0) {
        aplus1[nde, ] <- rnorm(nnd * nsim, mean = 0, sd = 1)
    }
    P1pd <- array(P1[nde,nde],c(nnd,nnd))
nde<-array(nde,c(nnd))
storage.mode(nnd) <- "integer"
storage.mode(nde) <-"integer"
storage.mode(P1pd) <- "double"
    sims.out <- .Fortran("simsmoother", PACKAGE = "KFAS", NAOK = TRUE, 
        ydimt, tv, y, Z, Tt, Rt, H, Qt, a1, P1, P1pd=P1pd, P1inf, nnd, nde, nsim, 
        alphasim = alphasim, epsplus, etaplus, aplus1, p, 
        n, m, r, info = info, tol)
    if (sims.out$info != 0) {
        if (sims.out$info == 1) 
            stop("Couldn't compute Cholesky factorization of Ht!")
        if (sims.out$info == 2) 
            stop("Couldn't compute Cholesky factorization of Qt!")
        if (sims.out$info == 3) 
            stop("Couldn't compute Cholesky factorization of P1!")
        else stop("Error in filtering!")
    }
    return(sims.out$alphasim)
}
