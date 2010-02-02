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
   
    p <- dim(yt)[1]
    n <- dim(yt)[2]
    m <- length(a1)
    if (is.vector(Qt)) 
        r <- 1
    else r <- dim(as.array(Qt))[2]
    tv <- array(0, dim = 5)
    tv[1] <- !(is.na(dim(as.array(Tt))[3]) || dim(as.array(Tt))[3] == 
        1)
    tv[2] <- !(is.na(dim(as.array(Rt))[3]) || dim(as.array(Rt))[3] == 
        1)
    tv[3] <- !(is.na(dim(as.array(Qt))[3]) || dim(as.array(Qt))[3] == 
        1)
tv[4] <- !(is.na(dim(as.array(Ht))[3]) || dim(as.array(Ht))[3] == 
        1)
tv[5] <- !(is.na(dim(as.array(Zt))[3]) || dim(as.array(Zt))[3] == 
        1)

    ymiss <- is.na(yt)
    ydimt <- array(p, dim = n)

tv[5]<-max(tv[4],tv[5])

if(sum(ymiss)>0)
{
	tv[4]<-tv[5]<-1
}
    H <- array(Ht, c(p, p, (n-1)*tv[4]+1))
    Z <- array(Zt, dim = c(p, m, (n-1)*tv[5]+1))
    y <- yt
if(sum(ymiss)>0){
    for (i in 1:n) {
        ydimt[i] <- sum(!ymiss[1:p, i])
        if (ydimt[i] != p && ydimt[i] != 0) {
            y[1:ydimt[i], i] <- yt[!ymiss[, i], i]
            H[1:ydimt[i], 1:ydimt[i], i] <- H[!ymiss[, i], !ymiss[, 
                i], i]
            Z[1:ydimt[i], , i] <- Z[!ymiss[, i], , i]
        }
    }}
    at <- array(0, dim = c(m, n + 1))
    Pt <- Pinf <- Pstar <- array(0, dim = c(m, m, n + 1))    
    Kt <- Ktuni <- Kinf <- Kstar <- Kinfuni <- Kstaruni <- array(0, dim = c(m, p, n))     
    vt <- vtuni <- Ftuni <- Finfuni <- Fstaruni <- array(0, dim = c(p, n))  
    Ft <- Finf <- Fstar <- array(0, dim = c(p, p, n))    
    Lt <- Linf <- Lstar <- array(0, dim = c(m, m, n))
     

    Pinf[, , 1] <- P1inf
    lik <- info <- j <- d <- 0
    storage.mode(d) <- storage.mode(j) <- storage.mode(p) <- storage.mode(m) <- storage.mode(r) <- storage.mode(n) <- storage.mode(tv) <- storage.mode(info) <- storage.mode(optcal) <- storage.mode(ydimt) <- storage.mode(j) <- "integer"

    kfout <- NULL  
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

kfout$tv[4] <- !(is.na(dim(as.array(Ht))[3]) || dim(as.array(Ht))[3] == 
        1)
kfout$tv[5]  <- !(is.na(dim(as.array(Zt))[3]) || dim(as.array(Zt))[3] == 
        1)
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

    kfout$yt <- yt
    kfout$Tt <- Tt
    kfout$Rt <- Rt
    kfout$Qt <- Qt
    kfout$Zt <- Zt
    kfout$Ht <- Ht
Zt<-array(Zt,c(p,m,(n-1)*kfout$tv[5]+1))
Ht<-array(Ht,c(p,p,(n-1)*kfout$tv[4]+1))
  if (kfout$d > 0) {
	kfout$Kt[, , 1:kfout$d] <- NA
        kfout$Ft[, , 1:kfout$d] <- NA
        kfout$Pt[, , 1:kfout$d] <- NA
        for (i in 1:kfout$d) {
            if (ydimt[i] != p) {
		kfout$Finfuni[(ydimt[i] + 1):p, i] <- kfout$Fstaruni[(ydimt[i] + 1):p, i] <- kfout$vtuni[(ydimt[i] + 1):p, i] <- NA
		if(sum(optcal)>0)
		{
			if(optcal[2]==1){
				kfout$Fstar[, , i] <- matrix(Zt[, , (i-1)*kfout$tv[5]+1],p,m) %*% kfout$Pstar[, , i] %*% t(matrix(Zt[, , (i-1)*kfout$tv[5]+1],p,m)) + Ht[, , (i-1)*kfout$tv[4]+1]
        		        kfout$Finf[, , i] <- matrix(Zt[, , (i-1)*kfout$tv[5]+1],p,m) %*% kfout$Pinf[, , i] %*% t(matrix(Zt[,,(i-1)*kfout$tv[5]+1],p,m))
			}                
        	        if (optcal[1]) {
        			kfout$vt[!ymiss[,i],i]<-kfout$vt[1:ydimt[i], i]
				kfout$vt[ymiss[,i],i]<-NA}
			if (optcal[3]) {
				kfout$Kinf[,!is.na(yt[,i]), i]<-kfout$Kinf[,1:kfout$ydimt[i], i]
				kfout$Kinf[,is.na(yt[,i]), i]<-NA
				kfout$Kstar[,!is.na(yt[,i]), i]<-kfout$Kstar[,1:kfout$ydimt[i], i]
				kfout$Kstar[,is.na(yt[,i]), i]<-NA
			}
            	}
	  }
    }}
    if (kfout$d < n) {
        for (i in (kfout$d + 1):n) {
            if (ydimt[i] != p) {
		kfout$Ftuni[(ydimt[i] + 1):p, i] <- kfout$vtuni[(ydimt[i] + 1):p, i] <- NA  
		if(sum(optcal)>0){
			if(optcal[2]==1){
				kfout$Ft[, , i] <- matrix(Zt[, ,(i-1)*kfout$tv[5]+1],p,m) %*% kfout$Pt[, , i] %*% t(matrix(Zt[,,(i-1)*kfout$tv[5]+1],p,m)) + Ht[, , (i-1)*kfout$tv[4]+1]
			}          
                if (optcal[1]) {
		  kfout$vt[!ymiss[,i],i]<-kfout$vt[1:ydimt[i], i]
		  kfout$vt[ymiss[,i],i]<-NA}    
		if (optcal[3]) {
				kfout$Kt[,!is.na(yt[,i]), i]<-kfout$Kt[,1:kfout$ydimt[i], i]
				kfout$Kt[,is.na(yt[,i]), i]<-NA				
			}        
            }}
        }
    }
 if (optcal[1] == 0) 
        kfout$vt <- NULL
    if (optcal[2] == 0) {
        kfout$Ft <- kfout$Finf <- kfout$Fstar <- NULL
    }
    else{
	kfout$Finf <- array(kfout$Finf[, , 1:kfout$d], c(p, p, kfout$d))
    	kfout$Fstar <- array(kfout$Fstar[, , 1:kfout$d], c(p, p, kfout$d))	
	}
    if (optcal[3] == 0) {
        kfout$Kt <- kfout$Kinf <- kfout$Kstar <- NULL
    }
    else{
	kfout$Kinf <- array(kfout$Kinf[, , 1:kfout$d], c(m, p, kfout$d))
    	kfout$Kstar <- array(kfout$Kstar[, , 1:kfout$d], c(m, p, kfout$d))
	}
    if (optcal[4] == 0) {
        kfout$Lt <- kfout$Linf <- kfout$Lstar <- NULL
    }
    else{
    	kfout$Linf <- array(kfout$Linf[, , 1:kfout$d], c(m, m, kfout$d))
    	kfout$Lstar <- array(kfout$Lstar[, , 1:kfout$d], c(m, m, kfout$d))
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
    Nt0 <- Nt1 <- Nt2 <- array(0, dim = c(out$m, out$m, out$d + 1))
    rt <- array(0, dim = c(out$m, out$n + 1))
    rt0 <-rt1 <- array(0, dim = c(out$m, out$d + 1))
    ymiss <- is.na(out$yt)
tv<-out$tv
tv[5]<-max(tv[4],tv[5])
if(sum(ymiss)>0)
{
tv[4]<-tv[5]<-1
}
    H <- array(out$Ht, c(out$p, out$p, (out$n-1)*tv[4]+1))
    Z <- array(out$Zt, dim = c(out$p, out$m, (out$n-1)*tv[5]+1))
if(sum(ymiss)>0){
    for (i in 1:out$n) {
	 if (out$ydimt[i] != out$p && out$ydimt[i] != 0) {
            H[1:out$ydimt[i], 1:out$ydimt[i], i] <- H[!ymiss[, i], !ymiss[, 
                i], i]
            Z[1:out$ydimt[i], , i] <- Z[!ymiss[, i], , i]
        }
    }}

     ks.out <- .Fortran("ks", PACKAGE = "KFAS", NAOK = TRUE, out$ydimt, 
        tv, Z = Z, array(out$Tt, c(out$m, out$m, (out$n - 
            1) * out$tv[1] + 1)), H = H, out$at, 
        out$Pt, out$vtuni, out$Ftuni, out$Ktuni, ahat = ahat, 
        Vt = Vt, rt = rt, rt0 = rt0, rt1 = rt1, Nt = Nt, Nt0 = Nt0, 
        Nt1 = Nt1, Nt2 = Nt2, Pinf = out$Pinf, 
        Pstar = out$Pstar, Kinfuni = out$Kinfuni, Kstaruni = out$Kstaruni, 
        Finfuni = out$Finfuni, Fstaruni = out$Fstaruni, d = out$d, 
        j = out$j, p = out$p, m = out$m, n = out$n, 
        tol = out$tol)
    ks.out <- c(out, list(ahat = ks.out$ahat, Vt = ks.out$Vt, 
        rt = ks.out$rt, rt0 = ks.out$rt0, rt1 = ks.out$rt1, Nt = ks.out$Nt, 
        Nt0 = ks.out$Nt0, Nt1 = ks.out$Nt1, Nt2 = ks.out$Nt2)) 
       
    ks.out
}

distsmoother<-function(out)
{
epshat <- array(0, dim = c(out$p, out$n))
epshatvar <- array(0, dim = c(out$p, out$p, out$n))
etahat <- array(0, dim = c(out$r, out$n))
etahatvar <- array(0, dim = c(out$r, out$r, out$n))

tv <- array(0, dim = 3)
    tv[1] <- !(is.na(dim(as.array(out$Ht))[3]) || dim(as.array(out$Ht))[3] == 
        1)
    tv[2] <- !(is.na(dim(as.array(out$Rt))[3]) || dim(as.array(out$Rt))[3] == 
        1)
    tv[3] <- !(is.na(dim(as.array(out$Qt))[3]) || dim(as.array(out$Qt))[3] == 
        1)
storage.mode(tv)<-"integer"
info<-0
storage.mode(info)<-"integer"
ds.out<-.Fortran("distsmooth",PACKAGE="KFAS",NAOK=TRUE,tv, array(out$Ht, c(out$p, out$p, (out$n - 
            1) * tv[1] + 1)), array(out$Rt, c(out$m, out$r, (out$n - 1) * tv[2] + 1)), array(out$Qt, 
            c(out$r, out$r, (out$n - 1) * tv[3] + 1)),out$Ft,out$Kt,out$vt,out$Nt,out$rt,out$Fstar,out$Finf,out$Kinf,out$Kstar,out$Nt0,out$rt0,out$d,out$p,out$m,out$r,out$n,out$tol,epshat=epshat,epshatvar=epshatvar,etahat=etahat,etahatvar=etahatvar,info)
if(info==1) stop("Fstar singular!")
if(info==2) stop("Ft singular!")
c(out,list(epshat = ds.out$epshat, epshatvar = ds.out$epshatvar,etahat = ds.out$etahat, etahatvar = ds.out$etahatvar))
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
    storage.mode(yt) <-"double"
    p <- dim(yt)[1]
    n <- dim(yt)[2]
    m <- length(a1)
    if (is.vector(Qt)) 
        r <- 1
    else r <- dim(as.array(Qt))[2]
    tv <- array(0, dim = 5)
    tv[1] <- !(is.na(dim(as.array(Tt))[3]) || dim(as.array(Tt))[3] == 
        1)
    tv[2] <- !(is.na(dim(as.array(Rt))[3]) || dim(as.array(Rt))[3] == 
        1)
    tv[3] <- !(is.na(dim(as.array(Qt))[3]) || dim(as.array(Qt))[3] == 
        1)
 tv[4] <- !(is.na(dim(as.array(Ht))[3]) || dim(as.array(Ht))[3] == 
        1)
    tv[5] <- !(is.na(dim(as.array(Zt))[3]) || dim(as.array(Zt))[3] == 
        1)
    ymiss <- is.na(yt)
    ydimt <- array(p, dim = n)

tv[5]<-max(tv[4],tv[5])

if(sum(ymiss)>0)
{
	tv[4]<-tv[5]<-1
}
    H <- array(Ht, c(p, p, (n-1)*tv[4]+1))
    Z <- array(Zt, dim = c(p, m, (n-1)*tv[5]+1))
    y <- yt
if(sum(ymiss)>0){
    for (i in 1:n) {
        ydimt[i] <- sum(!ymiss[1:p, i])
        if (ydimt[i] != p && ydimt[i] != 0) {
            y[1:ydimt[i], i] <- yt[!ymiss[, i], i]
            H[1:ydimt[i], 1:ydimt[i], i] <- H[!ymiss[, i], !ymiss[, 
                i], i]
            Z[1:ydimt[i], , i] <- Z[!ymiss[, i], , i]
        }
    }}

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
    storage.mode(tv) <- storage.mode(ydimt) <- storage.mode(p) <- storage.mode(m) <- storage.mode(r) <- storage.mode(n) <- storage.mode(info) <- storage.mode(nsim) <- "integer"
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
    sims.out <- .Fortran("simsmoother", PACKAGE = "KFAS", NAOK = TRUE, 
        ydimt, tv, y, Z, Tt, Rt, H, Qt, a1, P1, P1pd, P1inf, nnd, nde, nsim, 
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
    return(alphasim=sims.out$alphasim)
}

eflik0<-function(yt,Zt,Tt,Rt,Qt,a1,P1,P1inf, dist=c("Poisson", "Binomial", "Negative binomial"), offset=1)
{                                          
                                   
dist <- match.arg(dist)
if(dist == "Poisson")
	distr<-1
else {
	if(dist == "Binomial")
		distr<-2
	else 
		distr<-3
}
storage.mode(distr)<-"integer"
n<-length(yt)
p<-1
m <- length(a1)
if (is.vector(Qt)) 
        r <- 1
else r <- dim(as.array(Qt))[2]

alpha<-matrix(0,m,n)
tv <- array(0, dim = 3)
    tv[1] <- !(is.na(dim(as.array(Tt))[3]) || dim(as.array(Tt))[3] == 
        1)
    tv[2] <- !(is.na(dim(as.array(Rt))[3]) || dim(as.array(Rt))[3] == 
        1)
    tv[3] <- !(is.na(dim(as.array(Qt))[3]) || dim(as.array(Qt))[3] == 
        1)
 tv[4] <- 1
    tv[5] <- !(is.na(dim(as.array(Zt))[3]) || dim(as.array(Zt))[3] == 
        1)
   
       
ydimt <- array(c(!is.na(yt)), dim = n)

if(sum(ydimt) < n)
{
tv[5]<-1
}

at <- rt <- rt0 <- rt1 <- array(0, dim = c(m, n + 1))
vt <- vtuni <- Ftuni <- Finfuni <- Fstaruni <- epshat <- ytilde <-array(0, dim = c(1, n))
Ft <- Finf <- Fstar <- epshatvar <- Ht <- array(0, dim = c(1, 1, n))
Pt <- Pinf <- Pstar <- Nt <- Nt0 <- Nt1 <- Nt2 <- array(0, dim = c(m, m, n + 1))
Kt <- Ktuni <- Kinf <- Kstar <- Kinfuni <- Kstaruni <- array(0, dim = c(m, 1, n))
Lt <- Linf <- Lstar <- Vt <- array(0, dim = c(m, m, n))
Pinf[, , 1] <- P1inf
lik <- info <- d <- j <- 0
ahat <- array(0, dim = c(m, n))
etahat <- array(0, dim = c(r, n))
etahatvar <- array(0, dim = c(r, r, n))
theta <- array(0,dim=n)
offset <- array(offset,dim=n)
optcal<-array(1,dim=4)
storage.mode(distr)<-storage.mode(n)<-storage.mode(p)<-storage.mode(r)<-storage.mode(m)<-storage.mode(d)<-storage.mode(j)<-storage.mode(tv)<-storage.mode(ydimt)<-storage.mode(info)<-storage.mode(optcal)<-"integer"

out<-.Fortran("eflik0", PACKAGE = "KFAS", NAOK = TRUE,  yt = array(yt,dim=c(1,n)), 
        ydimt = ydimt, tv = tv, Zt = array(Zt,c(1,m,(n-1)*tv[5]+1)), Tt = array(Tt,c(m,m,(n-1)*tv[1]+1)), Rt = array(Rt,c(m,r,(n-1)*tv[2]+1)), Ht = Ht, Qt = array(Qt,c(r,r,(n-1)*tv[3]+1)), a1 = a1, P1 = P1, at = at, Pt = Pt, vtuni = vtuni, Ftuni = Ftuni, 
        Ktuni = Ktuni, Pinf = Pinf, Pstar = Pstar, Finfuni = Finfuni, 
        Fstaruni = Fstaruni, Kinfuni = Kinfuni, Kstaruni = Kstaruni, 
        d = d, j = j, p = p, m = m, r = r, n = n, lik = lik, 
        optcal = optcal, info = info, vt = vt, Ft = Ft, Kt = Kt, 
        Lt = Lt, Finf = Finf, Fstar = Fstar, Kinf = Kinf, Kstar = Kstar, 
        Linf = Linf, Lstar = Lstar, ahat = ahat, Vt = Vt, rt = rt, rt0 = rt0, rt1 = rt1, Nt = Nt, Nt0 = Nt0, 
        Nt1 = Nt1, Nt2 = Nt2, epshat=epshat, epshatvar=epshatvar, etahat=etahat, etahatvar=etahatvar, tol = 1e-7, theta=theta, offset=offset,ytilde=ytilde,dist=distr)

if(dist=="Poisson")
{
	ueth<- offset*exp(out$theta)
	lik <- out$lik + sum(dpois(out$yt[1,],ueth,log=TRUE) - dnorm(out$ytilde[1,],out$theta,sqrt(out$Ht[1,1,]),log=TRUE)) + log(1 - .125*sum(ueth*out$epshatvar[1,1,]^2))
}else {
	if(dist=="Binomial")
		{
			ueth<- exp(out$theta)/(1+exp(out$theta))
			lik <- out$lik + sum(dbinom(x=out$yt[1,],size=offset,prob=ueth,log=TRUE) - dnorm(out$ytilde[1,],out$theta,sqrt(out$Ht[1,1,]),log=TRUE)) + log(1 - .125*sum(offset*(1-ueth^2)/(ueth-2)^3))
		}else
			{
				ueth <- 1 - exp(out$theta)
				lik <- out$lik + sum(dnbinom(x=out$yt[1,],size=offset,prob=ueth,log=TRUE) - dnorm(out$ytilde[1,],out$theta,sqrt(out$Ht[1,1,]),log=TRUE)) + log(1 - .125*sum(offset*exp(out$theta)*(exp(2*out$theta)+4*exp(out$theta)+1)/(exp(out$theta)-4)^4))
			}
}
out$P1inf<-P1inf
out$lik0<-lik
out$dist<-dist
invisible(out)
}


efsmoother <-function(out,nsim) {

  alphasim <-simsmoother(out$ytilde, out$Zt, out$Tt, out$Rt, out$Ht, out$Qt, out$a1, out$P1, out$P1inf, nsim)
  n<-out$n
  m<-out$m
  thetasim<-array(0,dim=c(1,n,nsim))

Zt<-array(out$Zt,c(1,m,n))

for (k in 1:nsim) {
    for (t in 1:n) {
      thetasim[,t,k] <- Zt[,,t]%*%alphasim[,t,k]
    }
  }  


  alphahat <- array(0,dim=c(m,n))
  alphavar <- array(0,dim=c(m,m,n))

  thetahat <- array(0,dim=c(1,n))
  thetavar <- array(0,dim=c(1,1,n))

  division <- vector(length=nsim)

if(out$dist=="Poisson")
{
for(k in 1:nsim)
{
#division[k] <- prod(dpois(out$yt[1, ], out$offset * exp(thetasim[1, , k])))/prod(dnorm(out$ytilde[1, ], mean = thetasim[1, , k], sd = sqrt(out$Ht[1,1, ]))) 
division[k] <-exp(sum(dpois(out$yt[1,],out$offset*exp(thetasim[1, , k]),log=TRUE)-dnorm(out$ytilde[1,],mean=thetasim[1, , k],sd=sqrt(out$Ht[1,1,]),log=TRUE)))

}
}
else{ if(out$dist=="Binomial")
{
for(k in 1:nsim)
{
division[k] <- prod(dbinom(x=out$yt[1,],size=offset,prob=exp(thetasim[1,,k])/(1+exp(thetasim[1,,k]))))/prod(dnorm(out$ytilde[1,],mean=thetasim[1,,k],sd=sqrt(out$Ht[1,1,])))
}
}
else
for(k in 1:nsim)
{
division[k] <- prod(dnbinom(x=out$yt[1,],size=offset,prob=1 - exp(thetasim[1,,k])))/prod(dnorm(out$ytilde[1,],mean=thetasim[1,,k],sd=sqrt(out$Ht[1,1,])))
}
}


  denominator <- sum(division)

  for (t in 1:n) {

    numerator1 <- 0
    numerator2 <- 0
  
    for (k in 1:nsim) {
      numerator1 <- numerator1 + alphasim[,t,k]*division[k]
      numerator2 <- numerator2 + (alphasim[,t,k])%*%t(alphasim[,t,k])*division[k]
    }

    alphahat[,t] <- numerator1/denominator
    alphavar[,,t] <- numerator2/denominator - alphahat[,t]%*%t(alphahat[,t])
  }
 

  for (t in 1:n) {
    thetahat[,t] <- Zt[,,t]%*%alphahat[,t]
  }

  for (t in 1:n) {
    thetavar[,,t] <- Zt[,,t]%*%alphavar[,,t]%*%t(t(Zt[,,t]))
  }  
 
  
  
  out <- c(out,list(nsim=nsim, ahat=alphahat,
              that=thetahat,
              ahatvar=alphavar,
              thatvar=thetavar))
  invisible(out)

}







eflik <-function(yt,Zt,Tt,Rt,Qt,a1,P1,P1inf, dist=c("Poisson", "Binomial", "Negative binomial"), offset=1,nsim=1000){

dist <- match.arg(dist)
if(dist == "Poisson")
	distr<-1
else {
	if(dist == "Binomial")
		distr<-2
	else 
		distr<-3
}
storage.mode(distr)<-"integer"
n<-length(yt)
p<-1
m <- length(a1)
if (is.vector(Qt)) 
        r <- 1
else r <- dim(as.array(Qt))[2]

alpha<-matrix(0,m,n)
tv <- array(0, dim = 5)
    tv[1] <- !(is.na(dim(as.array(Tt))[3]) || dim(as.array(Tt))[3] == 
        1)
    tv[2] <- !(is.na(dim(as.array(Rt))[3]) || dim(as.array(Rt))[3] == 
        1)
    tv[3] <- !(is.na(dim(as.array(Qt))[3]) || dim(as.array(Qt))[3] == 
        1)
tv[4] <- 1
    tv[5] <- !(is.na(dim(as.array(Zt))[3]) || dim(as.array(Zt))[3] == 
        1)

ydimt <- array(c(!is.na(yt)), dim = n)

if(sum(ydimt) < n)
{
tv[5]<-1
}

at <- rt <- rt0 <- rt1 <- array(0, dim = c(m, n + 1))
vt <- vtuni <- Ftuni <- Finfuni <- Fstaruni <- ytilde <-array(0, dim = c(1, n))
Ft <- Finf <- Fstar <- Ht <- array(0, dim = c(1, 1, n))
Pt <- Pinf <- Pstar <- Nt <- Nt0 <- Nt1 <- Nt2 <- array(0, dim = c(m, m, n + 1))
Kt <- Ktuni <- Kinf <- Kstar <- Kinfuni <- Kstaruni <- array(0, dim = c(m, 1, n))
Lt <- Linf <- Lstar <- Vt <- array(0, dim = c(m, m, n))
Pinf[, , 1] <- P1inf
lik <- info <- d <- j <- 0
ahat <- array(0, dim = c(m, n))
theta <- array(0,dim=n)
offset <- array(offset,dim=n)
optcal<-array(1,dim=4)
storage.mode(distr)<-storage.mode(n)<-storage.mode(p)<-storage.mode(r)<-storage.mode(m)<-storage.mode(d)<-storage.mode(j)<-storage.mode(tv)<-storage.mode(ydimt)<-storage.mode(info)<-storage.mode(optcal)<-"integer"


ytilde <- array(0, dim = c(p, n))

out<-.Fortran("eflik", PACKAGE = "KFAS", NAOK = TRUE,  yt = array(yt,dim=c(1,n)), 
        ydimt = ydimt, tv = tv, Zt = array(Zt,c(1,m,(n-1)*tv[5]+1)), Tt = array(Tt,c(m,m,(n-1)*tv[1]+1)), Rt = array(Rt,c(m,r,(n-1)*tv[2]+1)), Ht = Ht, Qt = array(Qt,c(r,r,(n-1)*tv[3]+1)), a1 = a1, P1 = P1, at = at, Pt = Pt, vtuni = vtuni, Ftuni = Ftuni, 
        Ktuni = Ktuni, Pinf = Pinf, Pstar = Pstar, Finfuni = Finfuni, 
        Fstaruni = Fstaruni, Kinfuni = Kinfuni, Kstaruni = Kstaruni, 
        d = d, j = j, p = p, m = m, r = r, n = n, lik = lik, 
        optcal = optcal, info = info, vt = vt, Ft = Ft, Kt = Kt, 
        Lt = Lt, Finf = Finf, Fstar = Fstar, Kinf = Kinf, Kstar = Kstar, 
        Linf = Linf, Lstar = Lstar, ahat = ahat, Vt = Vt, rt = rt, rt0 = rt0, rt1 = rt1, Nt = Nt, Nt0 = Nt0, 
        Nt1 = Nt1, Nt2 = Nt2, tol = 1e-7, theta=theta, offset=offset,ytilde=ytilde,dist=distr)


alphasim <-simsmoother(out$ytilde, out$Zt, out$Tt, out$Rt, out$Ht, out$Qt, a1, P1, P1inf, nsim)

thetasim<-array(0,dim=c(1,n,nsim))

Zt<-array(out$Zt,c(1,m,n))

for (k in 1:nsim) {
    for (t in 1:n) {
      thetasim[,t,k] <- Zt[,,t]%*%alphasim[,t,k]
    }
  }  

eg <- vector(length=nsim)

if(dist=="Poisson")
{
for(k in 1:nsim)
	{
	eg[k] <-sum(dpois(out$yt[1,],out$offset*exp(thetasim[1, , k]))/dnorm(out$ytilde[1,],mean=thetasim[1, , k],sd=sqrt(out$Ht[1,1,])))
	}
}
else
{ 
	if(dist=="Binomial")
	{
		for(k in 1:nsim)
		{
			eg[k] <- sum(dbinom(x=out$yt[1,],size=offset,prob=exp(thetasim[1,,k])/(1+exp(thetasim[1,,k])))/dnorm(out$ytilde[1,],mean=thetasim[1,,k],sd=sqrt(out$Ht[1,1,])))
		}
	}
	else
		for(k in 1:nsim)
		{
			eg[k] <- sum(dnbinom(x=out$yt[1,],size=offset,prob=1 - exp(thetasim[1,,k]))/dnorm(out$ytilde[1,],mean=thetasim[1,,k],sd=sqrt(out$Ht[1,1,])))
		}
}
out$likp<-out$lik+log(sum(eg))
out$dist<-dist
invisible(out)

invisible(out)

}


