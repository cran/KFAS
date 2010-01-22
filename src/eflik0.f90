subroutine eflik0(yt, ydimt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1, at, pt, vtuni, ftuni,& 
     ktuni, pinf, pstar, finfuni, fstaruni, kinfuni, kstaruni, d, j, p, m, r, n, lik,&
     optcal, info, vt, ft, kt, lt, finf, fstar, kinf, kstar, linf, lstar, ahat, vvt, &
     rt, rt0, rt1, nt, nt0, nt1, nt2, epshat, epshatvar, etahat, etahatvar, eps, theta, offset, ytilde,dist)

implicit none


integer ::  yna, tvh, tvhz
integer, intent(inout), dimension(1,n) :: ymiss

integer, intent(in) ::  p, m, r, n,dist
integer, intent(inout) :: d, j, info
integer :: i
double precision, intent(in) :: eps
integer, intent(inout), dimension(n) :: ydimt
integer, intent(inout), dimension(5) :: timevar
double precision, intent(in), dimension(1,n) :: yt
double precision, intent(in), dimension(1,m,n) :: zt  
double precision, intent(in), dimension(m,m,(n-1)*timevar(1)+1) :: tt 
double precision, intent(in), dimension(m,r,(n-1)*timevar(2)+1) :: rtv 
double precision, intent(inout), dimension(1,1,n) :: ht 
double precision, intent(in), dimension(r,r,(n-1)*timevar(3)+1) :: qt
double precision, intent(in), dimension(m) :: a1
double precision, intent(in), dimension(m,m) ::  p1
double precision, intent(inout), dimension(m,n+1) :: at,rt,rt0,rt1
double precision, intent(inout), dimension(1,n) :: vt,vtuni,ftuni,fstaruni,finfuni,epshat
double precision, intent(inout), dimension(1,1,n) :: ft,epshatvar,fstar,finf
double precision, intent(inout), dimension(m,m,n) :: lt,linf,lstar,vvt
double precision, intent(inout),dimension(m,m,n+1) ::  pt,pstar,pinf,nt,nt0,nt1,nt2
double precision, intent(inout),dimension(m,1,n) ::  kstaruni,kinfuni,kstar,kinf,kt,ktuni
double precision, intent(inout), dimension(m,n) :: ahat
double precision, intent(inout), dimension(r,n) :: etahat
double precision, intent(inout), dimension(r,r,n) :: etahatvar
double precision, intent(inout), dimension(n) :: theta,offset
double precision, intent(inout) :: lik
integer, intent(inout), dimension(4) :: optcal 
double precision :: err
double precision, dimension(m,n) :: alpha
double precision, dimension(1,n) :: ytilde
integer, dimension(3):: timevardis
double precision, dimension(1,n) ::  vtdis
double precision, dimension(1,1,n) :: ftdis
double precision, dimension(1,1,n) ::  fstardis

double precision, external :: ddot
external kf
external ks
external distsmooth

yna=0 
tvh=1
tvhz=1

err = 1.0d0
alpha = 0.0d0
optcal = 0 !

do while(err > 1d-7)
   do i=1,n
      theta(i) = ddot(m,zt(1,1:m,i),1,alpha(1:m,i),1)
   end do

   select case(dist)
   
      case(1)
         ht(1,1,1:n) = exp(-theta)/offset
         ytilde(1,1:n) = theta(1:n) + yt(1,1:n)*ht(1,1,1:n) - 1.0d0
         !ytilde(1,1:n) = theta(1:n) + ht(1,1,1:n)*(yt(1,1:n)-exp(theta))
      case(2)
         ht(1,1,1:n) = (1+exp(theta))**2/(offset*exp(theta))
         ytilde(1,1:n) = theta(1:n) + ht(1,1,1:n)*yt(1,1:n) - 1 - exp(theta)
      case(3)
         ht(1,1,1:n) = (1-exp(theta))**2/(offset*exp(theta))
         ytilde(1,1:n) = theta(1:n) + ht(1,1,1:n)*yt(1,1:n) + 1 - exp(-theta)
   end select

   at=0.0d0
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
   rt=0.0d0
   vvt=0.0d0
   rt0=0.0d0
   rt1=0.0d0
   nt0=0.0d0
   nt1=0.0d0
   nt2=0.0d0
   d=0
   j=0
ftdis = ht


call kf(ytilde, ymiss, ydimt, yna, tvh, tvhz, timevar, zt, tt, rtv, ftdis, qt, a1, p1, &
     at, pt, vtuni, ftuni, ktuni, pinf, pstar, finfuni, fstaruni, kinfuni, kstaruni, d, j, &
     p, m, r, n, lik, optcal, info, vt, ft, kt, lt, finf, fstar, kinf, kstar, linf, lstar, eps)

 ftdis = ht       


call ks(ymiss, yna,tvh,tvhz, timevar, zt, tt, ftdis, at, pt, vtuni, ftuni, ktuni, ahat, vvt, &
     rt, rt0, rt1, nt, nt0, nt1, nt2, pinf, pstar, kinfuni,&
     kstaruni, finfuni, fstaruni, d, j, p, m, n, eps)

!call ks(ymiss, yna,tvz,tvh,tvhz, timevar, zt, tt, ftdis, at, pt, vtuni, ftuni, ktuni, ahat, vvt, &
!     rt, rt0(1:m,1:(d+1)), rt1(1:m,1:(d+1)), nt, nt0(1:m,1:m,1:(d+1)), nt1(1:m,1:m,1:(d+1)), &
!     nt2(1:m,1:m,1:(d+1)), pinf(1:m,1:m,1:(d+1)), pstar(1:m,1:m,1:(d+1)), kinfuni(1:m,1:p,1:d),&
!     kstaruni(1:m,1:p,1:d), finfuni(1:p,1:d), fstaruni(1:p,1:d), d, j, p, m, n, eps)


   if(maxval(abs(ahat)) > 1d100) then
      info=10
      return
   end if
   err = maxval(abs(abs(alpha)-abs(ahat)))
   alpha = ahat
end do

optcal = 1

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
lik=0.0d0
d=0
j=0

ftdis = ht
lik=0.0d0


call kf(ytilde, ymiss, ydimt, yna, tvh, tvhz, timevar, zt, tt, rtv, ftdis, qt, a1, p1, &
     at, pt, vtuni, ftuni, ktuni, pinf, pstar, finfuni, fstaruni, kinfuni, kstaruni, d, j, &
     p, m, r, n, lik, optcal, info, vt, ft, kt, lt, finf, fstar, kinf, kstar, linf, lstar, eps)


ftdis = ht

call ks(ymiss, yna,tvh,tvhz, timevar, zt, tt, ftdis, at, pt, vtuni, ftuni, ktuni, ahat, vvt, &
     rt, rt0, rt1, nt, nt0, nt1, nt2, pinf, pstar, kinfuni,&
     kstaruni, finfuni, fstaruni, d, j, p, m, n, eps)

!call ks(ymiss, yna,tvz,tvh,tvhz, timevar, zt, tt, ftdis, at, pt, vtuni, ftuni, ktuni, ahat, vvt, &
!     rt, rt0(1:m,1:(d+1)), rt1(1:m,1:(d+1)), nt, nt0(1:m,1:m,1:(d+1)), nt1(1:m,1:m,1:(d+1)), &
!     nt2(1:m,1:m,1:(d+1)), pinf(1:m,1:m,1:(d+1)), pstar(1:m,1:m,1:(d+1)), kinfuni(1:m,1:p,1:d),&
!     kstaruni(1:m,1:p,1:d), finfuni(1:p,1:d), fstaruni(1:p,1:d), d, j, p, m, n, eps)

alpha = ahat
do i=1,n
   theta(i) = ddot(m,zt(1,1:m,i),1,alpha(1:m,i),1) !oli timevar
end do

timevardis(1) = 1
timevardis(2) = timevar(2)
timevardis(3) = timevar(3)
vtdis = vt
ftdis = ft
fstardis = fstar

call distsmooth(timevardis, ht, rtv, qt, ftdis, kt, vtdis, nt, rt, fstardis(1:p,1:p,1:d), finf(1:p,1:p,1:d),&
     kinf(1:m,1:p,1:d), kstar(1:m,1:p,1:d), nt0(1:m,1:m,1:(d+1)), rt0(1:m,1:(d+1)), d, p, m, r, n,&
     eps, epshat, epshatvar, etahat, etahatvar, info)



end subroutine
