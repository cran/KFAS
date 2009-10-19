subroutine simsmoother(ydimt, timevar, yt, zt, tt, rtv, ht, qt, a1, p1, p1pd, p1inf, &
     nnd, nde, nsim, alphasim, epsplus, etaplus, aplus1, p, n, m, r, info, eps)

  implicit none
  
  integer, intent(inout) :: info
  integer, intent(in) :: nsim, p, m, r, n, nnd
  integer ::  t, i, d, j
  integer, intent(inout), dimension(nnd) :: nde
  integer, intent(inout), dimension(n) :: ydimt
  integer, intent(in), dimension(3) :: timevar
  double precision, intent(inout), dimension(p,n) :: yt
  double precision, intent(inout), dimension(p,m,n) :: zt  
  double precision, intent(in), dimension(m,m,(n-1)*timevar(1)+1) :: tt 
  double precision, intent(in), dimension(m,r,(n-1)*timevar(2)+1) :: rtv 
  double precision, intent(inout), dimension(p,p,n) :: ht 
  double precision, intent(in), dimension(r,r,(n-1)*timevar(3)+1) :: qt
  double precision, intent(in), dimension(m) :: a1
  double precision, intent(in), dimension(m,m) ::  p1
  double precision, intent(in), dimension(nnd,nnd) ::  p1pd
  double precision, intent(in), dimension(m,m) ::  p1inf
  double precision, dimension(p,m,n) :: zthelp
  double precision, dimension(p,p,n) :: hthelp
  double precision, dimension(m,n+1) :: at
  double precision, dimension(m,m,n+1) :: pt
  double precision, dimension(p,n) :: vtuni
  double precision, dimension(p,n) :: vt
  double precision, dimension(p,p,n) :: ft
  double precision, dimension(p,n) :: ftuni
  double precision, dimension(m,p,n) :: kt
  double precision, dimension(m,p,n) :: ktuni
  double precision, dimension(m,m,n) :: lt
  double precision, dimension(m,m,n+1) ::  pstar
  double precision, dimension(m,m,n+1) ::  pinf
  double precision, dimension(m,p,n) ::  kstaruni
  double precision, dimension(m,p,n) ::  kinfuni
  double precision, dimension(m,p,n) ::  kstar
  double precision, dimension(m,p,n) ::  kinf
  double precision, dimension(p,n) ::  fstaruni
  double precision, dimension(p,n) ::  finfuni
  double precision, dimension(p,p,n) ::  fstar
  double precision, dimension(p,p,n) ::  finf
  double precision, dimension(m,m,n) :: linf
  double precision, dimension(m,m,n) :: lstar
  double precision :: lik
  integer, dimension(4) :: optcal
  double precision, dimension(m,m,n+1) :: nt !n_1 = n_0, ..., n_201 = n_200
  double precision, dimension(m,n+1) :: rt !same as n, r_1 = r_0 etc.
  double precision, dimension(m,n) :: ahat
  double precision, dimension(m,m,n) :: vvt
  double precision, dimension(m,n+1) :: rt0
  double precision, dimension(m,n+1) :: rt1
  double precision, dimension(m,m,n+1) :: nt0
  double precision, dimension(m,m,n+1) :: nt1
  double precision, dimension(m,m,n+1) :: nt2
  double precision, intent(inout), dimension(p,n,nsim) :: epsplus
  double precision, intent(inout), dimension(r,n,nsim) :: etaplus
  double precision, dimension(m,n+1) :: aplus
  double precision, intent(inout), dimension(m,nsim) :: aplus1
  double precision, dimension(m,n) :: aplushat
  double precision, dimension(p,n) :: yplus
  double precision, intent(inout), dimension(m,n,nsim) :: alphasim
  double precision, intent(in) :: eps
  double precision, dimension(p,p) :: cholht
  double precision, dimension(r,r) :: cholqt
  double precision, dimension(nnd,nnd) :: cholp1
  double precision, dimension(nnd) :: a1nd

  external dcopy
  external dgemm
  external dgemv
  external daxpy
  external dsyr
  external dger
  external dsymm
  external dsymv
  external dsyevr
  external dtrsm
  external dtrsv
  external kf
  external ks

  optcal = 0
  pinf=0.0d0 
  pinf(1:m,1:m,1) = p1inf(1:m,1:m)
  info = 0

  at=0.d0
  pt=0.0d0
  vtuni=0.0d0
  vt=0.0d0
  ft=0.0d0
  ftuni=0.0d0
  kt=0.0d0
  ktuni=0.0d0
  lt=0.0d0
  pstar=0.0d0
  kstaruni=0.0d0
  kinfuni=0.0d0
  kstar=0.0d0
  kinf=0.0d0
  fstaruni=0.0d0
  finfuni=0.0d0
  fstar=0.0d0
  finf=0.0d0
  linf=0.0d0
  lstar=0.0d0
  lik=0.0d0 
  nt=0.0d0
  rt =0.0d0
  ahat=0.0d0
  vvt=0.0d0
  rt0=0.0d0
  rt1=0.0d0
  nt0=0.0d0
  nt1=0.0d0
  nt2=0.0d0
  aplus=0.0d0  
  aplushat=0.0d0
  yplus=0.0d0
  alphasim=0.0d0
  d=0
  j=0
  
  zthelp = zt
  hthelp = ht


 call kf(yt, ydimt, timevar, zthelp, tt, rtv, hthelp, qt, a1, p1, at, pt, vtuni,&
          ftuni, ktuni, pinf, pstar, finfuni, fstaruni, kinfuni, kstaruni, d, j, p, m, r, n,&
          lik, optcal, info, vt, ft, kt, lt, finf, fstar, kinf, kstar, linf, lstar, eps)
 if(info /= 0) then
    info=4
    return
 end if
     
 zthelp = zt
 hthelp = ht

 call ks(ydimt, timevar, zthelp, tt, hthelp, at, pt, vtuni, ftuni, ktuni, ahat, vvt, &
          rt, rt0, rt1, nt, nt0, nt1, nt2, pinf, pstar, kinfuni,&
          kstaruni, finfuni, fstaruni, d, j, p, m, n, eps)

 do i = 1, floor(nsim/2.0d0)

    zthelp = zt
    hthelp = ht     
     
    do t = 1, n       
        if(ydimt(t)>0) then
           cholht(1:ydimt(t),1:ydimt(t)) = ht(1:ydimt(t),1:ydimt(t),t)
           call dpotrf('l',ydimt(t),cholht(1:ydimt(t),1:ydimt(t)),ydimt(t),info)
           if(info /= 0) then
              info=1
              return
           end if
           call dtrmv('l','n','n',ydimt(t),cholht(1:ydimt(t),1:ydimt(t)),ydimt(t),epsplus(1:ydimt(t),t,i),1)
        end if

        cholqt(1:r,1:r) = qt(1:r,1:r,(t-1)*timevar(3)+1)
        call dpotrf('l',r,cholqt,r,info)
        if(info /= 0) then
           info=2
           return
        end if
        call dtrmv('l','n','n',r,cholqt,r,etaplus(1:r,t,i),1)    
     end do

     if(nnd>0) then        
        cholp1 = p1pd
        call dpotrf('l',nnd,cholp1,nnd,info)
        if(info /= 0) then
           info=3
           return
        end if
        a1nd = aplus1(nde,i)
        call dtrmv('l','n','n',nnd,cholp1,nnd,a1nd,1)
        aplus(1:m,1) = a1(1:m)
        aplus(nde,1) = aplus(nde,1) + a1nd 
     end if
     
     do t = 1, n
        if(ydimt(t)>0) then
           call dgemv('n',ydimt(t),m,1.0d0,zt(1:ydimt(t),1:m,t),ydimt(t),aplus(1:m,t),1,0.0d0,yplus(1:ydimt(t),t),1)
           yplus(1:ydimt(t),t) = yplus(1:ydimt(t),t) + epsplus(1:ydimt(t),t,i)
        end if
        call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(1)+1),m,aplus(1:m,t),1,0.0d0,aplus(1:m,t+1),1)
        call dgemv('n',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,etaplus(1:r,t,i),1,1.0d0,aplus(1:m,t+1),1)
     end do
    
     pinf=0.0d0 
     pinf(1:m,1:m,1) = p1inf(1:m,1:m)
     info = 0
     
     at=0.d0
     pt=0.0d0
     vtuni=0.0d0
     ftuni=0.0d0
     ktuni=0.0d0
     pstar=0.0d0
     kstaruni=0.0d0
     kinfuni=0.0d0
     fstaruni=0.0d0
     finfuni=0.0d0
     nt=0.0d0
     rt =0.0d0
     vvt=0.0d0
     rt0=0.0d0
     rt1=0.0d0
     nt0=0.0d0
     nt1=0.0d0
     nt2=0.0d0
     d=0
     j=0

     call kf(yplus, ydimt, timevar, zthelp, tt, rtv, hthelp, qt, a1, p1, at, pt, vtuni,&
          ftuni, ktuni, pinf, pstar, finfuni, fstaruni, kinfuni, kstaruni, d, j, p, m, r, n,&
          lik, optcal, info, vt, ft, kt, lt, finf, fstar, kinf, kstar, linf, lstar, eps)

     if(info /= 0) then
        info=4
        return
     end if

zthelp = zt
hthelp = ht

     call ks(ydimt, timevar, zthelp, tt, hthelp, at, pt, vtuni, ftuni, ktuni, aplushat, vvt, &
          rt, rt0, rt1, nt, nt0, nt1, nt2, pinf, pstar, kinfuni,&
          kstaruni, finfuni, fstaruni, d, j, p, m, n, eps)

     do t = 1, n
        alphasim(1:m,t,i) = ahat(1:m,t) - aplushat(1:m,t) + aplus(1:m,t)
        alphasim(1:m,t,i+floor(nsim/2.0d0)) = ahat(1:m,t) + aplushat(1:m,t) - aplus(1:m,t)
     end do
  end do

if(mod(nsim,2) /= 0) then

zthelp = zt
hthelp = ht

   do t = 1, n       
      if(ydimt(t)>0) then
         cholht(1:ydimt(t),1:ydimt(t))=ht(1:ydimt(t),1:ydimt(t),t)
         call dpotrf('l',ydimt(t),cholht(1:ydimt(t),1:ydimt(t)),ydimt(t),info)
         if(info /= 0) then
            info=1
            return
         end if
         call dtrmv('l','n','n',ydimt(t),cholht(1:ydimt(t),1:ydimt(t)),ydimt(t),epsplus(1:ydimt(t),t,nsim),1)
      end if
      
      cholqt=qt(1:r,1:r,(t-1)*timevar(3)+1)
      call dpotrf('l',r,cholqt,r,info)
      if(info /= 0) then
         info=2
         return
      end if
      call dtrmv('l','n','n',r,cholqt,r,etaplus(1:r,t,nsim),1)    
   end do
   
   if(nnd>0) then        
      cholp1 = p1pd
      call dpotrf('l',nnd,cholp1,nnd,info)
      if(info /= 0) then
         info=3
         return
      end if
      a1nd = aplus1(nde,nsim)
      call dtrmv('l','n','n',nnd,cholp1,nnd,a1nd,1)
      aplus(1:m,1) = a1(1:m)
      aplus(nde,1) = aplus(nde,1) + a1nd 
   end if
     
   
   do t = 1, n
      if(ydimt(t)>0) then
         call dgemv('n',ydimt(t),m,1.0d0,zt(1:ydimt(t),1:m,t),ydimt(t),aplus(1:m,t),1,0.0d0,yplus(1:ydimt(t),t),1)
         yplus(1:ydimt(t),t) = yplus(1:ydimt(t),t) + epsplus(1:ydimt(t),t,nsim)
      end if
      call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(1)+1),m,aplus(1:m,t),1,0.0d0,aplus(1:m,t+1),1)
      call dgemv('n',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,etaplus(1:r,t,nsim),1,1.0d0,aplus(1:m,t+1),1)
   end do
   
   pinf=0.0d0 
   pinf(1:m,1:m,1) = p1inf(1:m,1:m)
   info = 0
   
   at=0.d0
   pt=0.0d0
   vtuni=0.0d0
   ftuni=0.0d0
   ktuni=0.0d0
   pstar=0.0d0
   kstaruni=0.0d0
   kinfuni=0.0d0
   fstaruni=0.0d0
   finfuni=0.0d0
   nt=0.0d0
   rt =0.0d0
   vvt=0.0d0
   rt0=0.0d0
   rt1=0.0d0
   nt0=0.0d0
   nt1=0.0d0
   nt2=0.0d0
   d=0
   j=0
   
   call kf(yplus, ydimt, timevar, zthelp, tt, rtv, hthelp, qt, a1, p1, at, pt, vtuni,&
        ftuni, ktuni, pinf, pstar, finfuni, fstaruni, kinfuni, kstaruni, d, j, p, m, r, n,&
        lik, optcal, info, vt, ft, kt, lt, finf, fstar, kinf, kstar, linf, lstar, eps)
   if(info /= 0) then
      info=4
      return
   end if
   
zthelp = zt
hthelp = ht

   call ks(ydimt, timevar, zthelp, tt, hthelp, at, pt, vtuni, ftuni, ktuni, aplushat, vvt, &
        rt, rt0, rt1, nt, nt0, nt1, nt2, pinf, pstar, kinfuni,&
        kstaruni, finfuni, fstaruni, d, j, p, m, n, eps)

   do t = 1, n
      alphasim(1:m,t,nsim) = ahat(1:m,t) - aplushat(1:m,t) + aplus(1:m,t)
   end do
end if
   
end subroutine
