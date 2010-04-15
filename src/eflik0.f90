subroutine eflik0(yt, ydimt, timevar, zt, tt, rtv, ht, qt, a1, p1, at, pt, vtuni, ftuni,& 
     ktuni, pinf, pstar, finfuni, fstaruni, kinfuni, kstaruni, d, j, p, m, r, n, lik,&
     optcal, info, vt, ft, kt, lt, finf, fstar, kinf, kstar, linf, lstar, ahat, vvt, &
     rt, rt0, rt1, nt, nt0, nt1, nt2, epshat, epshatvar, etahat, etahatvar, eps, theta, &
     offset, ytilde, dist,maxiter)

implicit none

integer, intent(in) ::  p, m, r, n
integer, intent(inout) :: d, j,maxiter
integer, intent(inout), dimension(4) :: info
integer ::  t, i, k
integer, intent(in), dimension(n) :: ydimt
integer, intent(in), dimension(5) :: timevar
integer, intent(in) ::  dist
double precision, intent(inout), dimension(n) :: offset
double precision, intent(inout), dimension(n) :: theta
double precision, dimension(n) :: theta0
double precision, intent(in), dimension(1,n) :: yt
double precision, intent(inout), dimension(1,n) :: ytilde
double precision, intent(inout), dimension(1,m,(n-1)*timevar(5)+1) :: zt  
double precision, intent(in), dimension(m,m,(n-1)*timevar(1)+1) :: tt 
double precision, intent(in), dimension(m,r,(n-1)*timevar(2)+1) :: rtv 
double precision, intent(inout), dimension(1,1,(n-1)*timevar(4)+1) :: ht 
double precision, intent(in), dimension(r,r,(n-1)*timevar(3)+1) :: qt
double precision, intent(in), dimension(m) :: a1
double precision, intent(in), dimension(m,m) ::  p1
double precision, intent(inout), dimension(m,n+1) :: at
double precision, intent(inout), dimension(m,m,n+1) :: pt
double precision, intent(inout), dimension(1,n) :: vt
double precision, intent(inout), dimension(1,n) :: vtuni
double precision, intent(inout), dimension(1,1,n) :: ft
double precision, intent(inout), dimension(1,n) :: ftuni
double precision, intent(inout), dimension(m,1,n) :: kt
double precision, intent(inout), dimension(m,1,n) :: ktuni
double precision, intent(inout), dimension(m,m,n) :: lt
double precision, intent(inout),dimension(m,m,n+1) ::  pstar
double precision, intent(inout),dimension(m,m,n+1) ::  pinf
double precision, intent(inout),dimension(m,1,n) ::  kstaruni
double precision, intent(inout),dimension(m,1,n) ::  kinfuni
double precision, intent(inout),dimension(m,1,n) ::  kstar
double precision, intent(inout),dimension(m,1,n) ::  kinf
double precision, intent(inout),dimension(1,n) ::  fstaruni
double precision, intent(inout),dimension(1,n) ::  finfuni
double precision, intent(inout),dimension(1,1,n) ::  fstar
double precision, intent(inout),dimension(1,1,n) ::  finf
double precision, intent(inout), dimension(m,m,n) :: linf
double precision, intent(inout), dimension(m,m,n) :: lstar
double precision, intent(inout) :: lik
integer, intent(inout), dimension(4) :: optcal 
double precision, dimension(m) :: arec
double precision, dimension(m,m) :: prec
double precision, dimension(m,m) ::  psrec
double precision, dimension(m,m) ::  pirec
double precision, dimension(m,r) :: mr
double precision, dimension(p,m) :: pm
double precision, dimension(m,m) ::  im
double precision, dimension(m) ::  m1
double precision, dimension(m,m) ::  mm
double precision, dimension(m,1) ::  mp
double precision, intent(in) :: eps
double precision, intent(inout), dimension(1,n) :: epshat
double precision, intent(inout), dimension(1,1,n) :: epshatvar
double precision, intent(inout), dimension(r,n) :: etahat
double precision, intent(inout), dimension(r,r,n) :: etahatvar
double precision :: err
double precision, intent(inout), dimension(m,m,n+1) :: nt !n_1 = n_0, ..., n_201 = n_200
double precision, intent(inout), dimension(m,n+1) :: rt !same as n, r_1 = r_0 etc.
double precision, intent(inout), dimension(m,n) :: ahat
double precision, intent(inout), dimension(m,m,n) :: vvt
double precision, dimension(m,m) :: linfuni
double precision, dimension(m,m) :: l0
double precision, dimension(m,m) :: lstaruni
double precision, intent(inout), dimension(m,n+1) :: rt0
double precision, intent(inout), dimension(m,n+1) :: rt1
double precision, intent(inout), dimension(m,m,n+1) :: nt0
double precision, intent(inout), dimension(m,m,n+1) :: nt1
double precision, intent(inout), dimension(m,m,n+1) :: nt2
double precision, dimension(m,m) ::  lthelp
double precision, dimension(m,m) :: nrec
double precision, dimension(m,m) :: nrec1
double precision, dimension(m,m) :: nrec2
double precision, dimension(m) :: rrec
double precision, dimension(m) :: rrec1
double precision, dimension(m) :: rhelp
double precision, dimension(m,m) ::  mm2


integer, dimension(3):: timevardis
double precision, dimension(p,n) ::  vtdis
double precision, dimension(p,p,n) :: ftdis,fstardis


external dcopy
external dgemm
external dgemv
external daxpy
external dsyr
external dger
external dposv
external dtrsm
external dtrsv
external dsymm
external dsymv
external dsyevr


double precision, external :: ddot


err = 1.0d0
optcal = 0 !
theta0=0.0d0
k=0

do while(err > 1e-6 .AND. k < maxiter)
  k=k+1
 select case(dist)
   
      case(1)
         do i=1,n
            ht(1,1,i) = exp(-theta(i))/offset(i)
            ytilde(1,i) =  yt(1,i)*ht(1,1,i) + theta(i) - 1.0d0
         end do
      case(2)
         ht(1,1,1:n) = (1+exp(theta(1:n)))**2/(offset(1:n)*exp(theta(1:n)))
         ytilde(1,1:n) = theta(1:n) + ht(1,1,1:n)*yt(1,1:n) - 1 - exp(theta(1:n))
     case(3)
         ht(1,1,1:n) = (1-exp(theta(1:n)))**2/(offset(1:n)*exp(theta(1:n)))
         ytilde(1,1:n) = theta(1:n) + ht(1,1,1:n)*yt(1,1:n) + 1 - exp(-theta(1:n))
   end select


d=0
j=0
 info = 0
t=0
i=0
 at=0.0d0
pt=0.0d0
 vt=0.0d0
 vtuni=0.0d0
 ft=0.0d0
 ftuni=0.0d0
 kt=0.0d0
 ktuni=0.0d0
 lt=0.0d0
  pstar=0.0d0
 pinf(1:m,1:m,2:n)=0.0d0
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
optcal =0
arec =0.0d0
prec =0.0d0
psrec =0.0d0
pirec =0.0d0
mr =0.0d0
pm =0.0d0
m1 =0.0d0
mm =0.0d0
mp =0.0d0
epshat =0.0d0
epshatvar =0.0d0
etahat =0.0d0
etahatvar =0.0d0
nt  =0.0d0
rt =0.0d0
ahat =0.0d0
vvt =0.0d0
linfuni =0.0d0
l0 =0.0d0
lstaruni =0.0d0
rt0 =0.0d0
rt1 =0.0d0
nt0 =0.0d0
nt1 =0.0d0
nt2 =0.0d0
lthelp =0.0d0
nrec =0.0d0
nrec1 =0.0d0
nrec2 =0.0d0
rrec =0.0d0
rrec1 =0.0d0
rhelp =0.0d0
mm2 =0.0d0


im = 0.0d0
do i = 1, m
   im(i,i) = 1.0d0
end do

d=0

if(maxval(abs(pinf(1:m,1:m,1))) > eps) then

   pstar(1:m,1:m,1) = p1
   psrec = pstar(1:m,1:m,1)
   pirec = pinf(1:m,1:m,1)
   at(1:m,1) = a1
   arec = a1
   diffuse: do while(d < n)
      d = d+1
      do j=1, ydimt(d)
         call dsymv('u',m,1.0d0,psrec,m,zt(j,1:m,(d-1)*timevar(5)+1),1,0.0d0,m1,1) 
         fstaruni(j,d) = ddot(m,zt(j,1:m,(d-1)*timevar(5)+1),1,m1,1)  + ht(j,j,(d-1)*timevar(4)+1)
         call dsymv('u',m,1.0d0,pirec,m,zt(j,1:m,(d-1)*timevar(5)+1),1,0.0d0,m1,1)
         finfuni(j,d) = ddot(m,zt(j,1:m,(d-1)*timevar(5)+1),1,m1,1)! finf
         vtuni(j,d) = ytilde(j,d) - ddot(m,zt(j,1:m,(d-1)*timevar(5)+1),1,arec,1) !arec
         call dsymv('u',m,1.0d0,psrec,m,zt(j,1:m,(d-1)*timevar(5)+1),1,0.0d0,kstaruni(1:m,j,d),1) ! kstar_t,i = pstar_t,i*t(z_t,i)
         call dsymv('u',m,1.0d0,pirec,m,zt(j,1:m,(d-1)*timevar(5)+1),1,0.0d0,kinfuni(1:m,j,d),1) ! kinf_t,i = pinf_t,i*t(z_t,i)
         if (finfuni(j,d) > eps) then
            call daxpy(m,vtuni(j,d)/finfuni(j,d),kinfuni(1:m,j,d),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
            call dsyr('u',m,fstaruni(j,d)/(finfuni(j,d)**2),kinfuni(1:m,j,d),1,psrec,m) !psrec = psrec +  kinf*kinf'*fstar/finf^2
            call dsyr2('u',m,-1.0d0/finfuni(j,d),kstaruni(1:m,j,d),1,kinfuni(1:m,j,d),1,psrec,m) !psrec = psrec -(kstar*kinf'+kinf*kstar')/finf
          
            
            call dger(m,m,-1.0d0/finfuni(j,d),kinfuni(1:m,j,d),1,kinfuni(1:m,j,d),1,pirec,m)
            lik = lik - 0.5d0*log(finfuni(j,d))
         else
            call daxpy(m,vtuni(j,d)/fstaruni(j,d),kstaruni(1:m,j,d),1,arec,1) !a_rec = a_rec + kstar(:,i,t)*vt(:,t)/fstar(i,t)
            call dsyr('u',m,(-1.0d0)/fstaruni(j,d),kstaruni(1:m,j,d),1,psrec,m) !psrec = psrec -kstar*kstar'/fstar
            lik = lik - 0.5d0*(log(fstaruni(j,d)) + vtuni(j,d)**2/fstaruni(j,d))
         end if
         if(maxval(abs(pirec)) < eps) then
            exit diffuse
         end if
      end do
           
      call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(d-1)*timevar(1)+1),m,arec,1,0.0d0,at(1:m,d+1),1)  !at(:,t+1) = matmul(tt,a_rec)
      call dcopy(m,at(1:m,d+1),1,arec,1) ! a_rec = at(:,t+1)      
      call dsymm('r','u',m,m,1.0d0,psrec,m,tt(1:m,1:m,(d-1)*timevar(1)+1),m,0.0d0,mm,m)
      !call dgemm('n','n',m,m,m,1.0d0,tt(1:m,1:m,(d-1)*timevar(1)+1),m,psrec,m,0.0d0,mm,m)
      call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(d-1)*timevar(1)+1),m,0.0d0,pstar(1:m,1:m,d+1),m)
      if(m /= r) then
         if(r>1) then
            call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(d-1)*timevar(3)+1),r,rtv(1:m,1:r,(d-1)*timevar(2)+1),m,0.0d0,mr,m)
            call dgemm('n','t',m,m,r,1.0d0,mr,m,rtv(1:m,1:r,(d-1)*timevar(2)+1),m,1.0d0,pstar(1:m,1:m,d+1),m)   
         else
            call dger(m,m,qt(1,1,(d-1)*timevar(3)+1),rtv(1:m,1,(d-1)*timevar(2)+1),1,rtv(1:m,1,(d-1)*timevar(2)+1),1,&
                 pstar(1:m,1:m,d+1),m)
         end if
      else
         pstar(1:m,1:m,d+1) = pstar(1:m,1:m,d+1) + qt(1:r,1:r,(d-1)*timevar(3)+1)
      end if
      psrec = pstar(1:m,1:m,d+1) 
      call dsymm('r','u',m,m,1.0d0,pirec,m,tt(1:m,1:m,(d-1)*timevar(1)+1),m,0.0d0,mm,m) 
      !call dgemm('n','n',m,m,m,1.0d0,tt(1:m,1:m,(d-1)*timevar(1)+1),m,pirec,m,0.0d0,mm,m)
      call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(d-1)*timevar(1)+1),m,0.0d0,pinf(1:m,1:m,d+1),m)
      pirec = pinf(1:m,1:m,d+1)   
      j=ydimt(d) ! ????
if(maxval(abs(pirec)) < eps ) then !lisatty
         exit diffuse
      end if
   end do diffuse

!non-diffuse filtering begins
   prec = psrec
   do i = j+1, ydimt(d)     
      vtuni(i,d) = ytilde(i,d) - ddot(m,zt(i,1:m,(d-1)*timevar(5)+1),1,arec,1) !vtuni         
      call dsymv('u',m,1.0d0,prec,m,zt(i,1:m,(d-1)*timevar(5)+1),1,0.0d0,m1,1) ! p symmetric!
      ftuni(i,d) = ddot(m,zt(i,1:m,(d-1)*timevar(5)+1),1,m1,1) + ht(i,i,(d-1)*timevar(4)+1)
      if (abs(ftuni(i,d)) > eps) then !ftuni/=0
         call daxpy(m,(1.0d0)/ftuni(i,d),m1,1,ktuni(1:m,i,d),1) !ktuni
         call daxpy(m,vtuni(i,d),ktuni(1:m,i,d),1,arec,1) !a_rec = a_rec + ktuni(:,i,t)*vtuni(:,t)
         call dsyr('u',m,-ftuni(i,d),ktuni(1:m,i,d),1,prec,m) !p_rec = p_rec - ktuni*ktuni'*ftuni(i,t)        
      end if
      lik = lik - 0.5*(log(ftuni(i,d)) + vtuni(i,d)**2/ftuni(i,d))         
   end do
   
   call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(d-1)*timevar(1)+1),m,arec,1,0.0d0,at(1:m,d+1),1)  !at(:,t+1) = matmul(tt,a_rec)
   !pt(:,:,t+1) = matmul(matmul(tt,p_rec),transpose(tt)) + matmul(matmul(rt,qt),transpose(rt))
   call dsymm('r','u',m,m,1.0d0,prec,m,tt(1:m,1:m,(d-1)*timevar(1)+1),m,0.0d0,mm,m)
   call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(d-1)*timevar(1)+1),m,0.0d0,pt(1:m,1:m,d+1),m)
   if(m /= r) then
      if(r>1) then
         call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(d-1)*timevar(3)+1),r,rtv(1:m,1:r,(d-1)*timevar(2)+1),m,0.0d0,mr,m)
         call dgemm('n','t',m,m,r,1.0d0,mr,m,rtv(1:m,1:r,(d-1)*timevar(2)+1),m,1.0d0,pt(1:m,1:m,d+1),m)   
      else
         call dger(m,m,qt(1,1,(d-1)*timevar(3)+1),rtv(1:m,1,(d-1)*timevar(2)+1),1,rtv(1:m,1,(d-1)*timevar(2)+1),1,&
              pt(1:m,1:m,d+1),m)
      end if
   else
      pt(1:m,1:m,d+1) = pt(1:m,1:m,d+1) + qt(1:r,1:r,(d-1)*timevar(3)+1)
   end if
   
   call dcopy(m,at(1:m,d+1),1,arec,1) ! a_rec =at(:,t+1)
   prec = pt(1:m,1:m,d+1)
   pstar(1:m,1:m,d+1) = prec  
end if


!Non-diffuse filtering continues from t=d+1, i=1


if(d==0) then
   prec = p1
   arec = a1!   call dcopy(m,a1,1,arec,1)
   at(1:m,1) = a1 !call dcopy(m,a1,1,at(1:m,1),1) !at(:,1) = a1
   pt(1:m,1:m,1) = p1
end if
do t = d+1, n
   do i = 1, ydimt(t)     
      vtuni(i,t) = ytilde(i,t) - ddot(m,zt(i,1:m,(t-1)*timevar(5)+1),1,arec,1) !univariate vt           
      call dsymv('u',m,1.0d0,prec,m,zt(i,1:m,(t-1)*timevar(5)+1),1,0.0d0,m1,1) ! p symmetric!
      ftuni(i,t) = ddot(m,zt(i,1:m,(t-1)*timevar(5)+1),1,m1,1)  + ht(i,i,(t-1)*timevar(4)+1) !ftuni      
      if (abs(ftuni(i,t)) > eps) then !ft/=0
         call daxpy(m,(1.0d0)/ftuni(i,t),m1,1,ktuni(1:m,i,t),1) !kt kirja
         call daxpy(m,vtuni(i,t),ktuni(1:m,i,t),1,arec,1) !a_rec = a_rec + kt(:,i,t)*vt(:,t) 
         call dsyr('u',m,-ftuni(i,t),ktuni(1:m,i,t),1,prec,m) !p_rec = p_rec - kt*kt'*ft(i,i,t)
      end if
      lik = lik - 0.5*(log(ftuni(i,t)) + vtuni(i,t)**2/ftuni(i,t))
   end do
   
   call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(1)+1),m,arec,1,0.0d0,at(1:m,t+1),1)  !at(:,t+1) = matmul(tt,a_rec)
   
   call dsymm('r','u',m,m,1.0d0,prec,m,tt(1:m,1:m,(t-1)*timevar(1)+1),m,0.0d0,mm,m)
   call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-1)*timevar(1)+1),m,0.0d0,pt(1:m,1:m,t+1),m)
   if(m/=r) then
      if(r>1) then
         call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,0.0d0,mr,m)
         call dgemm('n','t',m,m,r,1.0d0,mr,m,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,1.0d0,pt(1:m,1:m,t+1),m)   
      else
         call dger(m,m,qt(1,1,(t-1)*timevar(3)+1),rtv(1:m,1,(t-1)*timevar(2)+1),1,rtv(1:m,1,(t-1)*timevar(2)+1),& 
              1,pt(1:m,1:m,t+1),m) 
      end if
   else
      pt(1:m,1:m,t+1) = pt(1:m,1:m,t+1) + qt(1:r,1:r,(t-1)*timevar(3)+1)
   end if
   call dcopy(m,at(1:m,t+1),1,arec,1) ! a_rec =at(:,t+1)
   prec = pt(1:m,1:m,t+1)
end do




im = 0.0d0
do i = 1, m
   im(i,i) = 1.0d0
end do

rrec = 0.0d0
nrec = 0.0d0
nt(1:m,1:m,n+1) = 0.0d0 !t goes from n+1 to 1, not from n to 0 !
rt(1:m,n+1) = 0.0d0

do t = n, d+1, -1 !do until diffuse starts
   do i = ydimt(t), 1 , -1       
      if(abs(ftuni(i,t)) > eps) then 
         lthelp = im
         call dger(m,m,-1.0d0,ktuni(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(5)+1),1,lthelp,m) !l = I -kz 
         call dgemv('t',m,m,1.0d0,lthelp,m,rrec,1,0.0d0,rhelp,1) 
         rrec = rhelp + vtuni(i,t)/ftuni(i,t)*zt(i,1:m,(t-1)*timevar(5)+1)    
         call dsymm('l','u',m,m,1.0d0,nrec,m,lthelp,m,0.0d0,mm,m) !n*l
         call dgemm('t','n',m,m,m,1.0d0,lthelp,m,mm,m,0.0d0,nrec,m) !n = l'nl
         call dger(m,m,(1.0d0/ftuni(i,t)),zt(i,1:m,(t-1)*timevar(5)+1),1,zt(i,1:m,(t-1)*timevar(5)+1),1,nrec,m) ! n = n+z'z/f  
      end if   
   end do

   call dcopy(m,rrec,1,rt(1:m,t),1) !r_t-1 = r_t,0
   nt(1:m,1:m,t) = nrec !n_t-1 = n_t,0
   call dcopy(m,at(1:m,t),1,ahat(1:m,t),1) !ahat = at
   call dsymv('u',m,1.0d0,pt(1:m,1:m,t),m,rt(1:m,t),1,1.0d0,ahat(1:m,t),1) !ahat = ahat+pt*r_t-1
   vvt(1:m,1:m,t) = pt(1:m,1:m,t)
   call dsymm('l','u',m,m,1.0d0,pt(1:m,1:m,t),m,nt(1:m,1:m,t),m,0.0d0,mm,m) !pt*n_t-1
   call dsymm('r','u',m,m,-1.0d0,pt(1:m,1:m,t),m,mm,m,1.0d0,vvt(1:m,1:m,t),m) !pt*n_t-1*pt
   if(t>1) then
      call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,rrec,1,0.0d0,rhelp,1) !r_t,p=t_t-1'*r_t+1
      rrec = rhelp
      call dsymm('l','u',m,m,1.0d0,nrec,m,tt(1:m,1:m,(t-2)*timevar(1)+1),m,0.0d0,mm,m) !n*t
      call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,mm,m,0.0d0,nrec,m) !n_t,p = t'nt
   end if
end do

if(d>0) then
   t=d 
rt0(1:m,d+1)=rt(1:m,d+1) 
   if(ydimt(t)>j) then
      do i = ydimt(t), (j+1) , -1  
         if(ftuni(i,t) > eps) then 
            lthelp = im
            call dger(m,m,-1.0d0,ktuni(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(5)+1),1,lthelp,m) !l = i -kz
            call dgemv('t',m,m,1.0d0,lthelp,m,rrec,1,0.0d0,rhelp,1) 
            rrec=rhelp
            call daxpy(m,vtuni(i,t)/ftuni(i,t),zt(i,1:m,(t-1)*timevar(5)+1),1,rrec,1)
            call dgemm('n','n',m,m,m,1.0d0,nrec,m,lthelp,m,0.0d0,mm,m) !n*l
            call dgemm('t','n',m,m,m,1.0d0,lthelp,m,mm,m,0.0d0,nrec,m) !n = l'nl
            call dger(m,m,(1.0d0)/ftuni(i,t),zt(i,1:m,(t-1)*timevar(5)+1),1,zt(i,1:m,(t-1)*timevar(5)+1),1,nrec,m) ! n = n+z'z/f      
         end if
      end do
      rrec1 = 0.0d0
      nrec1 = 0.0d0
      nrec2 = 0.0d0
   else
      rrec1 = 0.0d0
      nrec1 = 0.0d0
      nrec2 = 0.0d0
   end if  


   do i = j, 1, -1 
      if(finfuni(i,t)>eps) then
         linfuni = im            
         call dger(m,m,-1.0d0/finfuni(i,t),kinfuni(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(5)+1),1,linfuni,m) !linf
         rhelp = -kstaruni(1:m,i,t)
         call daxpy(m,fstaruni(i,t)/finfuni(i,t),kinfuni(1:m,i,t),1,rhelp,1)
         l0=0.0d0
         call dger(m,m,(1.0d0/finfuni(i,t)),rhelp,1,zt(i,1:m,(t-1)*timevar(5)+1),1,l0,m) !l0

         call dgemv('t',m,m,1.0d0,linfuni,m,rrec1,1,0.0d0,rhelp,1) !rt1
         call dcopy(m,rhelp,1,rrec1,1)
         call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
         call daxpy(m,(vtuni(i,t)/finfuni(i,t)),zt(i,1:m,(t-1)*timevar(5)+1),1,rrec1,1)

         call dgemv('t',m,m,1.0d0,linfuni,m,rrec,1,0.0d0,rhelp,1) !rt0 
         rrec = rhelp      

 
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
            call dgemm('n','n',m,m,m,1.0d0,mm,m,linfuni,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
            
            call dger(m,m,-1.0d0*fstaruni(i,t)/(finfuni(i,t)**2.0d0),zt(i,1:m,(t-1)*timevar(5)+1)&
                 ,1,zt(i,1:m,(t-1)*timevar(5)+1),1,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2
  
            call dsymm('l','u',m,m,1.0d0,nrec,m,l0,m,0.0d0,mm,m) !mm= nt0*l0
           
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2 + l0'*nt0*l0
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
            call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = nt2 + linf'*nt1*l0
            call dgemm('t','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1'*linf
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1'*linf hUOm ntrans

            call dgemm('n','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1*linf !!!!!!!!!!
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec1,m) !nt1 = linf'*mm
            call dger(m,m,(1.0d0)/finfuni(i,t),zt(i,1:m,(t-1)*timevar(5)+1),1,zt(i,1:m,(t-1)*timevar(5)+1),1,nrec1,m) 
            call dsymm('l','u',m,m,1.0d0,nrec,m,linfuni,m,0.0d0,mm,m) !mm= nt0*linf
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'nt1*linf + z'z/finf
            
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec,m) !nt0 = linf'*mm

        
      else
         lstaruni= im
         call dger(m,m,(-1.0d0)/fstaruni(i,t),kstaruni(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(5)+1),1,lstaruni,m) !lstar = I -Kstar*Z/Fstar         
         call dgemv('t',m,m,1.0d0,lstaruni,m,rrec,1,0.0d0,rhelp,1)
         rrec = rhelp
         call daxpy(m,vtuni(i,t)/fstaruni(i,t),zt(i,1:m,(t-1)*timevar(5)+1),1,rrec,1) !r0 = Z'vtuni/Fstar - Lstar'r0
         call dgemv('t',m,m,1.0d0,lstaruni,m,rrec1,1,0.0d0,rhelp,1)
         rrec1=rhelp         
         
         call dgemm('t','n',m,m,m,1.0d0,lstaruni,m,nrec,m,0.0d0,mm,m) !mm =lstar'*nt0
         call dgemm('n','n',m,m,m,1.0d0,mm,m,lstaruni,m,0.0d0,nrec,m) !nt0 = lstar'*nt0*lstar
         call dger(m,m,(1.0d0)/fstaruni(i,t),zt(i,1:m,(t-1)*timevar(5)+1),1,zt(i,1:m,(t-1)*timevar(5)+1),1,nrec,m)  !nt0 = z'z/fstar+lstar'*nt0*lstar
         call dgemm('n','n',m,m,m,1.0d0,nrec1,m,lstaruni,m,0.0d0,mm,m) !mm = nt1*lstar
         nrec1 = mm
         call dgemm('n','n',m,m,m,1.0d0,nrec2,m,lstaruni,m,0.0d0,mm,m) !mm = nt1*lstar
         nrec2 = mm
        
      end if
   end do

   rt0(1:m,t) = rrec
   rt1(1:m,t) = rrec1
   nt0(1:m,1:m,t) = nrec
   nt1(1:m,1:m,t) = nrec1
   nt2(1:m,1:m,t) = nrec2

   call dcopy(m,at(1:m,t),1,ahat(1:m,t),1) !ahat = at
   call dgemv('n',m,m,1.0d0,pstar(1:m,1:m,t),m,rt0(1:m,t),1,1.0d0,ahat(1:m,t),1) !ahat = at + pstar * rt0_t
   call dgemv('n',m,m,1.0d0,pinf(1:m,1:m,t),m,rt1(1:m,t),1,1.0d0,ahat(1:m,t),1) !ahat = at + pstar * rt0_t + pinf*rt1_t
     
   vvt(1:m,1:m,t) = pstar(1:m,1:m,t)
   call dgemm('n','n',m,m,m,1.0d0,pstar(1:m,1:m,t),m,nt0(1:m,1:m,t),m,0.0d0,mm,m) !mm = pstar*nt0
   call dgemm('n','n',m,m,m,-1.0d0,mm,m,pstar(1:m,1:m,t),m,1.0d0,vvt(1:m,1:m,t),m) !vvt = pstar - pstar*nt0*pstar 
   call dgemm('n','n',m,m,m,1.0d0,pinf(1:m,1:m,t),m,nt1(1:m,1:m,t),m,0.0d0,mm,m) !mm = pinf*nt1
   call dgemm('n','n',m,m,m,-1.0d0,mm,m,pstar(1:m,1:m,t),m,0.0d0,mm2,m) !mm2 = -pinf*nt1*pstar
   vvt(1:m,1:m,t) = vvt(1:m,1:m,t) + mm2 + transpose(mm2) !vvt = pstar - pstar*nt0*pstar  -pinf*nt1*pstar - t(pinf*nt1*pstar)
   call dgemm('n','n',m,m,m,1.0d0,pinf(1:m,1:m,t),m,nt2(1:m,1:m,t),m,0.0d0,mm,m) !mm = pinf*nt2
   call dgemm('n','n',m,m,m,-1.0d0,mm,m,pinf(1:m,1:m,t),m,1.0d0,vvt(1:m,1:m,t),m) !vvt = vvt - pinf*nt2*pinf

   call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,rrec,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
   rrec = rhelp
   call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,rrec1,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
   rrec1 = rhelp
   call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,nrec2,m,0.0d0,mm,m) !mm =t'*nt2
   call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(1)+1),m,0.0d0,nrec2,m) !nt2 = t'*nt2*t
   call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,nrec1,m,0.0d0,mm,m) !mm =t'*nt2
   call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(1)+1),m,0.0d0,nrec1,m) !nt2 = t'*nt2*t
   call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,nrec,m,0.0d0,mm,m) !mm =t'*nt2
   call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(1)+1),m,0.0d0,nrec,m) !nt2 = t'*nt2*t



   do t=(d-1), 1, -1
       do i = ydimt(t), 1, -1
         
         if(finfuni(i,t)> eps) then
            linfuni = im            
            call dger(m,m,-1.0d0/finfuni(i,t),kinfuni(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(5)+1),1,linfuni,m) !linf
            rhelp = -kstaruni(1:m,i,t)
            call daxpy(m,fstaruni(i,t)/finfuni(i,t),kinfuni(1:m,i,t),1,rhelp,1)
            l0=0.0d0
            call dger(m,m,(1.0d0/finfuni(i,t)),rhelp,1,zt(i,1:m,(t-1)*timevar(5)+1),1,l0,m) !l0

            call dgemv('t',m,m,1.0d0,linfuni,m,rrec1,1,0.0d0,rhelp,1) !rt1
            call dcopy(m,rhelp,1,rrec1,1)
            call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
            call daxpy(m,vtuni(i,t)/finfuni(i,t),zt(i,1:m,(t-1)*timevar(5)+1),1,rrec1,1)
            call dgemv('t',m,m,1.0d0,linfuni,m,rrec,1,0.0d0,rhelp,1) !rt0 
            rrec = rhelp
            
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
            call dgemm('n','n',m,m,m,1.0d0,mm,m,linfuni,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
            
            call dger(m,m,-1.0d0*fstaruni(i,t)/(finfuni(i,t)**2.0d0),&
                 zt(i,1:m,(t-1)*timevar(5)+1),1,zt(i,1:m,(t-1)*timevar(5)+1),1,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2
  
            call dsymm('l','u',m,m,1.0d0,nrec,m,l0,m,0.0d0,mm,m) !mm= nt0*l0
           
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2 + l0'*nt0*l0
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
            call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = nt2 + linf'*nt1*l0
            call dgemm('t','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1'*linf
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1'*linf hUOm ntrans
          
            call dgemm('n','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1*linf !!!!!!!!!!
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec1,m) !nt1 = linf'*mm
            call dger(m,m,(1.0d0)/finfuni(i,t),zt(i,1:m,(t-1)*timevar(5)+1),1,zt(i,1:m,(t-1)*timevar(5)+1),1,nrec1,m) 
            !nt1 = linf'nt1'linf + z'z/finf
            call dsymm('l','u',m,m,1.0d0,nrec,m,linfuni,m,0.0d0,mm,m) !mm= nt0*linf
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'nt1*linf + z'z/finf
            
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec,m) !nt0 = linf'*mm
            
           
         else
            lstaruni= im
            call dger(m,m,(-1.0d0)/fstaruni(i,t),kstaruni(1:m,i,t),1,zt(i,1:m,(t-1)*timevar(5)+1),1,lstaruni,m) !lstar = I -Kstar*Z/Fstar         
            call dgemv('t',m,m,1.0d0,lstaruni,m,rrec,1,0.0d0,rhelp,1) !oli beta 1.0d0!!!!... JA miinusmerkki
            rrec = rhelp
            call daxpy(m,vtuni(i,t)/fstaruni(i,t),zt(i,1:m,(t-1)*timevar(5)+1),1,rrec,1) !r0 = Z'vtuni/Fstar - Lstar'r0
            call dgemv('t',m,m,1.0d0,lstaruni,m,rrec1,1,0.0d0,rhelp,1)
            rrec1=rhelp           
            
            call dgemm('t','n',m,m,m,1.0d0,lstaruni,m,nrec,m,0.0d0,mm,m) !mm =lstar'*nt0
            call dgemm('n','n',m,m,m,1.0d0,mm,m,lstaruni,m,0.0d0,nrec,m) !nt0 = lstar'*nt0*lstar
            call dger(m,m,(1.0d0)/fstaruni(i,t),zt(i,1:m,(t-1)*timevar(5)+1),1,zt(i,1:m,(t-1)*timevar(5)+1),1,nrec,m)  !nt0 = z'z/fstar+lstar'*nt0*lstar
            call dgemm('n','n',m,m,m,1.0d0,nrec1,m,lstaruni,m,0.0d0,mm,m) !mm = nt1*lstar
            nrec1 = mm
            call dgemm('n','n',m,m,m,1.0d0,nrec2,m,lstaruni,m,0.0d0,mm,m) !mm = nt2*lstar
            nrec2 = mm         
            
         end if
      end do
          
      
      rt0(1:m,t) = rrec
      rt1(1:m,t) = rrec1
      nt0(1:m,1:m,t) = nrec
      nt1(1:m,1:m,t) = nrec1
      nt2(1:m,1:m,t) = nrec2
      
      call dcopy(m,at(1:m,t),1,ahat(1:m,t),1) !ahat = at
      call dgemv('n',m,m,1.0d0,pstar(1:m,1:m,t),m,rt0(1:m,t),1,1.0d0,ahat(1:m,t),1) !ahat = at + pstar * rt0_t
      call dgemv('n',m,m,1.0d0,pinf(1:m,1:m,t),m,rt1(1:m,t),1,1.0d0,ahat(1:m,t),1) !ahat = at + pstar * rt0_t + pinf*rt1_t 
      
      vvt(1:m,1:m,t) = pstar(1:m,1:m,t)
      call dgemm('n','n',m,m,m,1.0d0,pstar(1:m,1:m,t),m,nt0(1:m,1:m,t),m,0.0d0,mm,m) !mm = pstar*nt0
      call dgemm('n','n',m,m,m,-1.0d0,mm,m,pstar(1:m,1:m,t),m,1.0d0,vvt(1:m,1:m,t),m) !vvt = pstar - pstar*nt0*pstar 
      call dgemm('n','n',m,m,m,1.0d0,pinf(1:m,1:m,t),m,nt1(1:m,1:m,t),m,0.0d0,mm,m) !mm = pinf*nt1
      call dgemm('n','n',m,m,m,-1.0d0,mm,m,pstar(1:m,1:m,t),m,0.0d0,mm2,m) !mm2 = -pinf*nt1*pstar
      vvt(1:m,1:m,t) = vvt(1:m,1:m,t) + mm2 + transpose(mm2) !vvt = pstar - pstar*nt0*pstar  -pinf*nt1*pstar - t(pinf*nt1*pstar)
      call dgemm('n','n',m,m,m,1.0d0,pinf(1:m,1:m,t),m,nt2(1:m,1:m,t),m,0.0d0,mm,m) !mm = pinf*nt2
      call dgemm('n','n',m,m,m,-1.0d0,mm,m,pinf(1:m,1:m,t),m,1.0d0,vvt(1:m,1:m,t),m) !vvt = vvt - pinf*nt2*pinf
      
      if(t>1) then
         call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,rrec,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
         rrec = rhelp
         call dgemv('t',m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,rrec1,1,0.0d0,rhelp,1,1) !tarkiSta tOimivUUS!
         rrec1 = rhelp
         call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,nrec2,m,0.0d0,mm,m) !mm =t'*nt2
         call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(1)+1),m,0.0d0,nrec2,m) !nt2 = t'*nt2*t
         call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,nrec1,m,0.0d0,mm,m) !mm =t'*nt2
         call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(1)+1),m,0.0d0,nrec1,m) !nt2 = t'*nt2*t
         call dgemm('t','n',m,m,m,1.0d0,tt(1:m,1:m,(t-2)*timevar(1)+1),m,nrec,m,0.0d0,mm,m) !mm =t'*nt2
         call dgemm('n','n',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-2)*timevar(1)+1),m,0.0d0,nrec,m) !nt2 = t'*nt2*t
      end if

   end do
end if


  
 do i=1,n
      theta(i) = ddot(m,zt(1,1:m,(i-1)*timevar(5)+1),1,ahat(1:m,i),1)
   end do

err= maxval(abs(theta-theta0))

theta0=theta

end do

maxiter=k



optcal = 1



   do t=1, n
      if(ydimt(t)>0) then      
         vt(1,t)=vtuni(1,t)       
      end if
   end do
   do t= 1, d
      if(ydimt(t)>0) then         
         fstar(1,1,t) = fstaruni(1,t)
         finf(1,1,t) = finfuni(1,t)
      end if
   end do
   do t=d+1, n
      if(ydimt(t)>0) then       
         ft(1,1,t) = ftuni(1,t)
      end if
   end do


   mulkd:  do t=1,d
      if(ydimt(t)>0) then
         if(finf(1,1,t)<eps) then !finf=0
            call dgemm('n','t',m,1,m,1.0d0,pstar(1:m,1:m,t),m,zt(1,1:m,(t-1)*timevar(5)+1),p,0.0d0,&
                 mp(1:m,1),m) !mp = pstar*z'
            call dgemm('n','n',m,1,m,1.0d0/fstar(1,1,t),tt(1:m,1:m,(t-1)*timevar(1)+1),m,mp(1:m,1),m,0.0d0,&
                 kstar(1:m,1,t),m) !kstar = t*mp           
         else
            pm(1,1:m) = zt(1,1:m,(t-1)*timevar(5)+1)
           
            kinf(1:m,1,t) = zt(1,1:m,(t-1)*timevar(5)+1)/finf(1,1,t) !kinf = z'*inv(finf)
            
            call dgemm('n','n',m,1,m,1.0d0,pinf(1:m,1:m,t),m,kinf(1:m,1,t),m,0.0d0,&
                 mp(1:m,1),m) !mp = pinf*kinf*z'*inv(finf)
            call dgemm('n','n',m,ydimt(t),m,1.0d0,tt(1:m,1:m,(t-1)*timevar(1)+1),m,mp(1:m,1),&
                 m,0.0d0,kinf(1:m,1,t),m) !kinf = t*pinf*z'*inv(finf)
               
            call dgemm('n','t',m,1,m,1.0d0,pstar(1:m,1:m,t),m,zt(1,1:m,(t-1)*timevar(5)+1),1,0.0d0,&
                    mp(1:m,1),m) !mp = pstar*z'    
            call dgemm('n','n',m,1,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(1)+1),m,mp(1:m,1),m,&
                 0.0d0,kstar(1:m,1,t),m) !kstar = t*mp
            
            call dgemm('n','n',m,1,1,1.0d0,kinf(1:m,1,t),m,&
                 fstar(1,1,t),1,0.0d0,mp(1:m,1),m) !mp =kinf*fstar
            kstar(1:m,1,t) = (kstar(1:m,1,t) - mp(1:m,1))/finf(1,1,t) ! note the -, wrong in formula of Dk2003, corrected form in appendix!
            pm(1,1:m) = kstar(1:m,1,t)/finf(1,1,t)
           
         end if
      end if
   end do mulkd
   mulk: do t=d+1, n
      if(ydimt(t)>0) then
         call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(1)+1),m,ktuni(1:m,1,t),1,0.0d0,kt(1:m,1,t),1) !kt=t*ktuni
      end if
   end do mulk
   


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
