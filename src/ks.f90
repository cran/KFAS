subroutine ks(ydimt, timevar, zt, tt, ht, rtv, qt, at, pt, vtuni, ftuni, ktuni, ahat, vvt, &
     rt, rt0, rt1, nt, nt0, nt1, nt2, epshat, epshatvar, etahat, etahatvar, pinf, pstar, kinfuni,&
     kstaruni, finfuni, fstaruni, d, j, p, m, r, n, eps)

implicit none

integer, intent(in) :: d, j, p, m, r, n
integer :: t, i, k, l
integer, intent(in), dimension(n) :: ydimt
integer, intent(in), dimension(4) :: timevar
double precision, intent(inout), dimension(p,m,n) :: zt  
double precision, intent(in), dimension(m,m,(n-1)*timevar(1)+1) :: tt 
double precision, intent(inout), dimension(p,p,n) :: ht
double precision, intent(in), dimension(m,r,(n-1)*timevar(2)+1) :: rtv
double precision, intent(in), dimension(r,r,(n-1)*timevar(3)+1) :: qt
double precision, intent(in), dimension(m,n+1) :: at
double precision, intent(in), dimension(m,m,n+1) :: pt
double precision, intent(in), dimension(p,n) ::  vtuni
double precision, intent(in), dimension(p,n) :: ftuni
double precision, intent(in), dimension(m,p,n) :: ktuni
double precision, intent(inout), dimension(m,m,n+1) :: nt !n_1 = n_0, ..., n_201 = n_200
double precision, intent(inout), dimension(m,n+1) :: rt !same as n, r_1 = r_0 etc.
double precision, intent(inout), dimension(m,n) :: ahat
double precision, intent(inout), dimension(m,m,n) :: vvt
double precision, intent(inout), dimension(m,m,d+1) ::  pinf
double precision, intent(in), dimension(m,m,d+1) ::  pstar
double precision, intent(in),dimension(m,p,n) ::  kinfuni
double precision, intent(in),dimension(m,p,n) ::  kstaruni
double precision, intent(in),dimension(p,d) ::  fstaruni
double precision, intent(in), dimension(p,d) ::  finfuni
double precision, dimension(m,m) :: linfuni
double precision, dimension(m,m) :: l0
double precision, dimension(m,m) :: lstaruni
double precision, intent(inout), dimension(m,d+1) :: rt0
double precision, intent(inout), dimension(m,d+1) :: rt1
double precision, intent(inout), dimension(m,m,d+1) :: nt0
double precision, intent(inout), dimension(m,m,d+1) :: nt1
double precision, intent(inout), dimension(m,m,d+1) :: nt2
double precision, intent(inout), dimension(p,n) :: epshat
double precision, intent(inout), dimension(p,n) :: epshatvar
double precision, intent(inout), dimension(r,n) :: etahat
double precision, intent(inout), dimension(r,r,n) :: etahatvar
double precision, dimension(m,m) ::  lt
double precision, dimension(m,m) :: nrec
double precision, dimension(m,m) :: nrec1
double precision, dimension(m,m) :: nrec2
double precision, dimension(m) :: rrec
double precision, dimension(m) :: rrec1
double precision, dimension(m) :: rhelp
double precision, dimension(r) :: r1
double precision, dimension(m,r) :: mr
double precision, dimension(m,r) :: mr2
double precision, dimension(m,m) ::  im
double precision, dimension(m,m) ::  mm
double precision, dimension(m,m) ::  mm2
integer, dimension(n) :: hdiagtest
double precision, intent(in) :: eps
double precision, dimension(p) :: evalu
double precision, dimension(p,p) :: evect
double precision, dimension(26*p) :: work
double precision, dimension(10*p) :: iwork
integer, dimension(2*p) :: isuppz
integer :: ehelp, info
double precision, dimension(p,m) :: pm
external dcopy
external dger
external dgemv
external daxpy
external dsymm
external dsymv
external dgemm
external dtrsm
external dsyevr
double precision, external :: ddot

hdiagtest=0

if(p>1) then
   do t= 1, n
      test:  do i = 1, ydimt(t)
         do k = i+1, ydimt(t) 
            if(abs(ht(k,i,t)) > eps) then
               hdiagtest(t)=1
               exit test
            end if
         end do
      end do test
   end do
   if(sum(hdiagtest)/=0) then !LDL' decomposition, from Numerical Methods in Scientific computing, vol II, Dahlquist
 do t = 1, n
         if(hdiagtest(t)==1 .AND. ydimt(t)>0) then          
            info=0
            singtest: do i=1 , ydimt(t)
               if(ht(i,i,t)<=0) then
                  info=1
                  exit singtest
               end if
            end do singtest
            if(info==1) then
               call dsyevr('V','A','L',ydimt(t),ht(1:ydimt(t),1:ydimt(t),t),ydimt(t),0.0d0,0.0d0,0,0,0.0d0,&
                    ehelp,evalu(1:ydimt(t)),evect(1:ydimt(t),1:ydimt(t)),ydimt(t),isuppz(1:2*ydimt(t)),work,26*p,iwork,10*p,info)
               if(info/=0) then
                  info=1
                  return
               end if
               do i=1, ydimt(t)
                  ht(i,i,t) = evalu(i)
               end do
               hdiagtest(t)=2              
            else
               do k = 1, ydimt(t)
                  do i = k+1, ydimt(t)
                     ht(i,k,t) = ht(i,k,t)/ht(k,k,t)
                     do l = k+1, i, -1
                        ht(i,l,t) = ht(i,l,t) -  ht(k,k,t)*ht(i,k,t)*ht(l,k,t)
                     end do
                  end do
               end do
            end if
         end if
      end do

      do t = 1, n
         if(ydimt(t)>0 .AND. hdiagtest(t)/=0) then
            if(hdiagtest(t)==1) then               
               call dtrsm('l','l','n','u',ydimt(t),m,1.0d0,ht(1:ydimt(t),1:ydimt(t),t),ydimt(t),zt(1:ydimt(t),1:m,t),ydimt(t)) !solve z*=inv(L) * zt     
            else
               call dgemm('t','n',ydimt(t),m,ydimt(t),1.0d0,evect(1:ydimt(t),1:ydimt(t)),ydimt(t),zt(1:ydimt(t),1:m,t),&
                    ydimt(t),0.0d0,pm(1:ydimt(t),1:m),ydimt(t))
               zt(1:ydimt(t),1:m,t) = pm(1:ydimt(t),1:m)
            end if
         end if
      end do
   end if
end if


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
         lt = im
         call dger(m,m,-1.0d0,ktuni(1:m,i,t),1,zt(i,1:m,t),1,lt,m) !l = I -kz 
         call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1) 
         rrec = rhelp
         call daxpy(m,vtuni(i,t)/ftuni(i,t),zt(i,1:m,t),1,rrec,1)
         call dsymm('l','u',m,m,1.0d0,nrec,m,lt,m,0.0d0,mm,m) !n*l
         call dgemm('t','n',m,m,m,1.0d0,lt,m,mm,m,0.0d0,nrec,m) !n = l'nl
         call dger(m,m,(1.0d0/ftuni(i,t)),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec,m) ! n = n+z'z/f
         !call dsyr('u',m,(1.0d0/ftuni(i,t)),zt(i,1:m,(t-1)*timevar(1)+1),1,nrec,m) ! n = n+z'z/f      
         epshat(i,t) = ht(i,i,t)*(vtuni(i,t)/ftuni(i,t) - ddot(m,ktuni(1:m,i,t),1,rrec,1))
         call dgemv('n',m,m,1.0d0,nrec,m,ktuni(1:m,i,t),1,1.0d0,rhelp,1)
         epshatvar(i,t) = ht(i,i,t) - (ht(i,i,t)**2)*(1.0d0/ftuni(i,t) + ddot(m,ktuni(1:m,i,t),1,rhelp,1))  
      end if   
   end do

   call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,rrec,1,0.0d0,r1,1)
   call dsymv('u',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,r1,1,0.0d0,etahat(1:r,t),1)
   etahatvar(1:r,1:r,t) = qt(1:r,1:r,(t-1)*timevar(3)+1)
   call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,0.0d0,mr,m)
   call dgemm('n','n',m,r,m,1.0d0,nrec,m,mr,m,0.0d0,mr2,m)
   call dgemm('t','n',r,r,m,1.0d0,mr,m,mr2,m,1.0d0,etahatvar(1:r,1:r,t),r)   
 
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
   if(ydimt(t)>j) then
      do i = ydimt(t), (j+1) , -1  
         if(ftuni(i,t) > eps) then 
            lt = im
            call dger(m,m,-1.0d0,ktuni(1:m,i,t),1,zt(i,1:m,t),1,lt,m) !l = i -kz
            call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1) 
            rrec=rhelp
            call daxpy(m,vtuni(i,t)/ftuni(i,t),zt(i,1:m,t),1,rrec,1)
            ! call dsymm('l','u',m,m,1.0d0,nrec,m,lt,m,0.0d0,mm,m) !n*l
            call dgemm('n','n',m,m,m,1.0d0,nrec,m,lt,m,0.0d0,mm,m) !n*l
            call dgemm('t','n',m,m,m,1.0d0,lt,m,mm,m,0.0d0,nrec,m) !n = l'nl
            call dger(m,m,(1.0d0)/ftuni(i,t),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec,m) ! n = n+z'z/f      
            epshat(i,t) = ht(i,i,t)*(vtuni(i,t)/ftuni(i,t) - ddot(m,ktuni(1:m,i,t),1,rrec,1))
            call dgemv('n',m,m,1.0d0,nrec,m,ktuni(1:m,i,t),1,0.0d0,rhelp,1)
            epshatvar(i,t) = ht(i,i,t) - (ht(i,i,t)**2)*(1.0d0/ftuni(i,t) + ddot(m,ktuni(1:m,i,t),1,rhelp,1))
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
         call dger(m,m,-1.0d0/finfuni(i,t),kinfuni(1:m,i,t),1,zt(i,1:m,t),1,linfuni,m) !linf
         rhelp = -kstaruni(1:m,i,t)
         call daxpy(m,fstaruni(i,t)/finfuni(i,t),kinfuni(1:m,i,t),1,rhelp,1)
         l0=0.0d0
         call dger(m,m,(1.0d0/finfuni(i,t)),rhelp,1,zt(i,1:m,t),1,l0,m) !l0

         call dgemv('t',m,m,1.0d0,linfuni,m,rrec1,1,0.0d0,rhelp,1) !rt1
        ! rrec1 = rhelp
         call dcopy(m,rhelp,1,rrec1,1)
         call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
         call daxpy(m,(vtuni(i,t)/finfuni(i,t)),zt(i,1:m,t),1,rrec1,1)

         call dgemv('t',m,m,1.0d0,linfuni,m,rrec,1,0.0d0,rhelp,1) !rt0 
         rrec = rhelp      

         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
         call dgemm('n','n',m,m,m,1.0d0,mm,m,linfuni,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
         call dger(m,m,(-fstaruni(i,t)/(finfuni(i,t)**2.0d0)),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec2,m) 
         !nt2 = linf'nt2'linf + z'z*fstar/finf^2 !kOrJattu
         call dsymm('l','u',m,m,1.0d0,nrec,m,l0,m,0.0d0,mm,m) !mm= nt0*l0

         call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2 + l0'*nt0*l0
         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
         call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = nt2 + linf'*nt1*l0
         call dgemm('t','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1*linf
         call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1'*linf hUOm ntrans
      
         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec1,m) !nt1 = linf'*mm
         call dger(m,m,(1.0d0)/finfuni(i,t),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec1,m) 
!nt1 = linf'nt1'linf + z'z/finf
         call dsymm('l','u',m,m,1.0d0,nrec,m,linfuni,m,0.0d0,mm,m) !mm= nt0*linf
         !call dgemm('n','n',m,m,m,1.0d0,nrec,m,linfuni,m,0.0d0,mm,m) !mm= nt0*linf
         call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'nt1'linf + z'z/finf
      
         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec,m) !nt0 = linf'*mm

         epshat(i,t) = -ht(i,i,t)/finfuni(i,t)*ddot(m,kinfuni(1:m,i,t),1,rrec,1)
         call dgemv('n',m,m,1.0d0,nrec,m,kinfuni(1:m,i,t),1,0.0d0,rhelp,1)
         epshatvar(i,t) = ht(i,i,t) - (ht(i,i,t)**2)/(finfuni(i,t)**2)*ddot(m,kinfuni(1:m,i,t),1,rhelp,1)
      else
         lstaruni= im
         call dger(m,m,(-1.0d0)/fstaruni(i,t),kstaruni(1:m,i,t),1,zt(i,1:m,t),1,lstaruni,m) !lstar = I -Kstar*Z/Fstar         
         call dgemv('t',m,m,1.0d0,lstaruni,m,rrec,1,0.0d0,rhelp,1)
         rrec = rhelp
         call daxpy(m,vtuni(i,t)/fstaruni(i,t),zt(i,1:m,t),1,rrec,1) !r0 = Z'vtuni/Fstar - Lstar'r0
         call dgemv('t',m,m,1.0d0,lstaruni,m,rrec1,1,0.0d0,rhelp,1)
         rrec1=rhelp         
         
         call dgemm('t','n',m,m,m,1.0d0,lstaruni,m,nrec,m,0.0d0,mm,m) !mm =lstar'*nt0
         call dgemm('n','n',m,m,m,1.0d0,mm,m,lstaruni,m,0.0d0,nrec,m) !nt0 = lstar'*nt0*lstar
         call dger(m,m,(1.0d0)/fstaruni(i,t),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec,m)  !nt0 = z'z/fstar+lstar'*nt0*lstar
         call dgemm('n','n',m,m,m,1.0d0,nrec1,m,lstaruni,m,0.0d0,mm,m) !mm = nt1*lstar
         nrec1 = mm
         call dgemm('n','n',m,m,m,1.0d0,nrec2,m,lstaruni,m,0.0d0,mm,m) !mm = nt1*lstar
         nrec2 = mm
      
         epshat(i,t) = ht(i,i,t) * (vtuni(i,t)/fstaruni(i,t) - ddot(m,kinfuni(1:m,i,t),1,rrec,1)/fstaruni(i,t))
         call dgemv('n',m,m,1.0d0,nrec,m,kstaruni(1:m,i,t),1,0.0d0,rhelp,1)
         epshatvar(i,t) = ht(i,i,t) - ht(i,i,t)**2 * (1.0d0/fstaruni(i,t) + &
              ddot(m,kstaruni(1:m,i,t),1,rhelp,1)/(fstaruni(i,t)**2))
      end if
   end do
   call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,rrec,1,0.0d0,r1,1)
   call dsymv('u',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,r1,1,0.0d0,etahat(1:r,t),1)
   etahatvar(1:r,1:r,t) = qt(1:r,1:r,(t-1)*timevar(3)+1)
   call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,0.0d0,mr,m)
   call dgemm('n','n',m,r,m,1.0d0,nrec,m,mr,m,0.0d0,mr2,m)
   call dgemm('t','n',r,r,m,1.0d0,mr,m,mr2,m,1.0d0,etahatvar(1:r,1:r,t),r)

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
            call dger(m,m,-1.0d0/finfuni(i,t),kinfuni(1:m,i,t),1,zt(i,1:m,t),1,linfuni,m) !linf
            rhelp = -1.0d0*kstaruni(1:m,i,t)
            call daxpy(m,fstaruni(i,t)/finfuni(i,t),kinfuni(1:m,i,t),1,rhelp,1)
            l0=0.0d0
            call dger(m,m,1.0d0/finfuni(i,t),rhelp,1,zt(i,1:m,t),1,l0,m) !l0

            call dgemv('t',m,m,1.0d0,linfuni,m,rrec1,1,0.0d0,rhelp,1) !rt1
            rrec1 = rhelp
            call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
            call daxpy(m,vtuni(i,t)/finfuni(i,t),zt(i,1:m,t),1,rrec1,1)
            call dgemv('t',m,m,1.0d0,linfuni,m,rrec,1,0.0d0,rhelp,1) !rt0 
            rrec = rhelp
            
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
            call dgemm('n','n',m,m,m,1.0d0,mm,m,linfuni,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
            
            call dger(m,m,-1.0d0*fstaruni(i,t)/(finfuni(i,t)**2.0d0),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2
  
            call dsymm('l','u',m,m,1.0d0,nrec,m,l0,m,0.0d0,mm,m) !mm= nt0*l0
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2 + l0'*nt0*l0
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
            call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = nt2 + linf'*nt1*l0
            call dgemm('t','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1*linf
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1'*linf hUOm ntrans
            
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec1,m) !nt1 = linf'*mm
            call dger(m,m,(1.0d0)/finfuni(i,t),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec1,m) 
            !nt1 = linf'nt1'linf + z'z/finf
            call dsymm('l','u',m,m,1.0d0,nrec,m,linfuni,m,0.0d0,mm,m) !mm= nt0*linf
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'nt1'linf + z'z/finf
            
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec,m) !nt0 = linf'*mm
            
            epshat(i,t) = -ht(i,i,t)/finfuni(i,t)*ddot(m,kinfuni(1:m,i,t),1,rrec,1)
            call dgemv('n',m,m,1.0d0,nrec,m,kinfuni(1:m,i,t),1,1.0d0,rhelp,1)
            epshatvar(i,t) = ht(i,i,t) - (ht(i,i,t)**2)/(finfuni(i,t)**2)*&
                 (ftuni(i,t) + ddot(m,kinfuni(1:m,i,t),1,rhelp,1))
         else
            lstaruni= im
            call dger(m,m,(-1.0d0)/fstaruni(i,t),kstaruni(1:m,i,t),1,zt(i,1:m,t),1,lstaruni,m) !lstar = I -Kstar*Z/Fstar         
            call dgemv('t',m,m,1.0d0,lstaruni,m,rrec,1,0.0d0,rhelp,1) !oli beta 1.0d0!!!!... JA miinusmerkki
            rrec = rhelp
            call daxpy(m,vtuni(i,t)/fstaruni(i,t),zt(i,1:m,t),1,rrec,1) !r0 = Z'vtuni/Fstar - Lstar'r0
            call dgemv('t',m,m,1.0d0,lstaruni,m,rrec1,1,0.0d0,rhelp,1)
            rrec1=rhelp           
            
            call dgemm('t','n',m,m,m,1.0d0,lstaruni,m,nrec,m,0.0d0,mm,m) !mm =lstar'*nt0
            call dgemm('n','n',m,m,m,1.0d0,mm,m,lstaruni,m,0.0d0,nrec,m) !nt0 = lstar'*nt0*lstar
            call dger(m,m,(1.0d0)/fstaruni(i,t),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec,m)  !nt0 = z'z/fstar+lstar'*nt0*lstar
            call dgemm('n','n',m,m,m,1.0d0,nrec1,m,lstaruni,m,0.0d0,mm,m) !mm = nt1*lstar
            nrec1 = mm
            call dgemm('n','n',m,m,m,1.0d0,nrec2,m,lstaruni,m,0.0d0,mm,m) !mm = nt1*lstar
            nrec2 = mm
            
            epshat(i,t) = ht(i,i,t) * (vtuni(i,t)/fstaruni(i,t) - &
                 ddot(m,kinfuni(1:m,i,t),1,rrec,1)/fstaruni(i,t))
            call dgemv('n',m,m,1.0d0,nrec,m,kstaruni(1:m,i,t),1,0.0d0,rhelp,1)
            epshatvar(i,t) = ht(i,i,t) - ht(i,i,t)**2 * (1.0d0/fstaruni(i,t) &
                 + ddot(m,kstaruni(1:m,i,t),1,rhelp,1)/(fstaruni(i,t)**2))
         end if
      end do
      
      call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,rrec,1,0.0d0,r1,1)
      call dsymv('u',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,r1,1,0.0d0,etahat(1:r,t),1)
      etahatvar(1:r,1:r,t) = qt(1:r,1:r,(t-1)*timevar(3)+1)
      call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,0.0d0,mr,m)
      call dgemm('n','n',m,r,m,1.0d0,nrec,m,mr,m,0.0d0,mr2,m)
      call dgemm('t','n',r,r,m,1.0d0,mr,m,mr2,m,1.0d0,etahatvar(1:r,1:r,t),r)
      
      
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

end subroutine
