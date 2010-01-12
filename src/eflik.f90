subroutine eflik(yt, ydimt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1, at, pt, vtuni, ftuni,& 
     ktuni, pinf, pstar, finfuni, fstaruni, kinfuni, kstaruni, d, j, p, m, r, n, lik,&
     optcal, info, vt, ft, kt, lt, finf, fstar, kinf, kstar, linf, lstar, ahat, vvt, &
     rt, rt0, rt1, nt, nt0, nt1, nt2, eps, theta, offset, ytilde,dist)

implicit none



integer ::  yna, tvh, tvhz
integer, intent(in), dimension(1,n) :: ymiss
integer, intent(in) ::  p, m, r, n,dist
integer, intent(inout) :: d, j, info
integer :: i
double precision, intent(in) :: eps
integer, intent(in), dimension(n) :: ydimt
integer, intent(in), dimension(5) :: timevar
double precision, intent(in), dimension(p,n) :: yt
double precision, intent(in), dimension(p,m,(n-1)*timevar(5)+1) :: zt  
double precision, intent(in), dimension(m,m,(n-1)*timevar(1)+1) :: tt 
double precision, intent(in), dimension(m,r,(n-1)*timevar(2)+1) :: rtv 
double precision, intent(inout), dimension(p,p,n) :: ht 
double precision, intent(in), dimension(r,r,(n-1)*timevar(3)+1) :: qt
double precision, intent(in), dimension(m) :: a1
double precision, intent(in), dimension(m,m) ::  p1
double precision, intent(inout), dimension(m,n+1) :: at
double precision, intent(inout), dimension(m,m,n+1) :: pt
double precision, intent(inout), dimension(p,n) :: vt
double precision, intent(inout), dimension(p,n) :: vtuni
double precision, intent(inout), dimension(p,p,n) :: ft
double precision, intent(inout), dimension(p,n) :: ftuni
double precision, intent(inout), dimension(m,p,n) :: kt
double precision, intent(inout), dimension(m,p,n) :: ktuni
double precision, intent(inout), dimension(m,m,n) :: lt
double precision, intent(inout),dimension(m,m,n+1) ::  pstar
double precision, intent(inout),dimension(m,m,n+1) ::  pinf
double precision, intent(inout),dimension(m,p,n) ::  kstaruni
double precision, intent(inout),dimension(m,p,n) ::  kinfuni
double precision, intent(inout),dimension(m,p,n) ::  kstar
double precision, intent(inout),dimension(m,p,n) ::  kinf
double precision, intent(inout),dimension(p,n) ::  fstaruni
double precision, intent(inout),dimension(p,n) ::  finfuni
double precision, intent(inout),dimension(p,p,n) ::  fstar
double precision, intent(inout),dimension(p,p,n) ::  finf
double precision, intent(inout), dimension(m,m,n) :: linf
double precision, intent(inout), dimension(m,m,n) :: lstar
double precision, intent(inout), dimension(m,m,n+1) :: nt !n_1 = n_0, ..., n_201 = n_200
double precision, intent(inout), dimension(m,n+1) :: rt !same as n, r_1 = r_0 etc.
double precision, intent(inout), dimension(m,n) :: ahat
double precision, intent(inout), dimension(m,m,n) :: vvt
double precision, intent(inout), dimension(m,n+1) :: rt0
double precision, intent(inout), dimension(m,n+1) :: rt1
double precision, intent(inout), dimension(m,m,n+1) :: nt0
double precision, intent(inout), dimension(m,m,n+1) :: nt1
double precision, intent(inout), dimension(m,m,n+1) :: nt2
double precision, intent(inout), dimension(n) :: theta
double precision, intent(inout), dimension(n) :: offset
double precision, intent(inout) :: lik
integer, intent(inout), dimension(4) :: optcal 
double precision :: err
double precision, dimension(m,n) :: alpha
double precision, dimension(p,n) :: ytilde
!double precision, dimension(n)  :: ueth

double precision, dimension(p,p,n) :: ftdis

double precision, external :: ddot
external kf
external ks
external distsmooth

yna=0 
tvh=1
tvhz=0
err = 1.0d0
alpha = 0.0d0
optcal = 0 !
do while(err > 1e-6)
   do i=1,n
      theta(i) = ddot(m,zt(1,1:m,(i-1)*timevar(5)+1),1,alpha(1:m,i),1)
   end do

   select case(dist)
   
      case(1)
         ht(1,1,1:n) = exp(-theta)/offset
         ytilde(1,1:n) = theta(1:n) + yt(1,1:n)*ht(1,1,1:n) - 1.0d0
      case(2)
         ht(1,1,1:n) = (1+exp(theta))**2/(offset*exp(theta))
         ytilde(1,1:n) = theta(1:n) + ht(1,1,1:n)*yt(1,1:n) - 1 - exp(theta)
      case(3)
         ht(1,1,1:n) = (1-exp(theta))**2/(offset*exp(theta))
         ytilde(1,1:n) = theta(1:n) + ht(1,1,1:n)*yt(1,1:n) + 1 - exp(-theta)
   end select

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
ftdis = ht

call kf(ytilde, ymiss, ydimt, yna, tvh, tvhz, timevar, zt, tt, rtv, ftdis, qt, a1, p1, &
     at, pt, vtuni, ftuni, ktuni, pinf, pstar, finfuni, fstaruni, kinfuni, kstaruni, d, j, &
     p, m, r, n, lik, optcal, info, vt, ft, kt, lt, finf, fstar, kinf, kstar, linf, lstar, eps)

 ftdis = ht       


call ks(ymiss, yna,tvh,tvhz, timevar, zt, tt, ftdis, at, pt, vtuni, ftuni, ktuni, ahat, vvt, &
     rt, rt0(1:m,1:(d+1)), rt1(1:m,1:(d+1)), nt, nt0(1:m,1:m,1:(d+1)), nt1(1:m,1:m,1:(d+1)), &
     nt2(1:m,1:m,1:(d+1)), pinf(1:m,1:m,1:(d+1)), pstar(1:m,1:m,1:(d+1)), kinfuni(1:m,1:p,1:d),&
     kstaruni(1:m,1:p,1:d), finfuni(1:p,1:d), fstaruni(1:p,1:d), d, j, p, m, n, eps)

   if(maxval(abs(ahat)) > 1.0d+100) then
      info=1
      return
   end if
   err = maxval(abs(alpha-ahat))
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
     rt, rt0(1:m,1:(d+1)), rt1(1:m,1:(d+1)), nt, nt0(1:m,1:m,1:(d+1)), nt1(1:m,1:m,1:(d+1)), &
     nt2(1:m,1:m,1:(d+1)), pinf(1:m,1:m,1:(d+1)), pstar(1:m,1:m,1:(d+1)), kinfuni(1:m,1:p,1:d),&
     kstaruni(1:m,1:p,1:d), finfuni(1:p,1:d), fstaruni(1:p,1:d), d, j, p, m, n, eps)

alpha = ahat
do i=1,n
   theta(i) = ddot(m,zt(1,1:m,((i-1)*timevar(5)+1)),1,alpha(1:m,i),1)
end do


end subroutine
