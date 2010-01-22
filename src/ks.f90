subroutine ks(ymiss, yna,tvh,tvhz, timevar, zt, tt, ht, at, pt, vtuni, ftuni, ktuni, ahat, vvt, &
     rt, rt0, rt1, nt, nt0, nt1, nt2, pinf, pstar, kinfuni,&
     kstaruni, finfuni, fstaruni, d, j, p, m, n, eps)

implicit none

integer, intent(in) :: d, j, p, m, n,yna,tvh,tvhz
integer :: t, i,k,info
integer, intent(in), dimension(5) :: timevar
integer, intent(in), dimension(p,n) :: ymiss
double precision, intent(inout), dimension(p,m,n) :: zt  
double precision, intent(in), dimension(m,m,(n-1)*timevar(1)+1) :: tt 
double precision, intent(inout), dimension(p,p,(n-1)*timevar(4)+1) :: ht 
double precision, intent(in), dimension(m,n+1) :: at
double precision, intent(in), dimension(m,m,n+1) :: pt
double precision, intent(in), dimension(p,n) ::  vtuni
double precision, intent(in), dimension(p,n) :: ftuni
double precision, intent(in), dimension(m,p,n) :: ktuni
double precision, intent(inout), dimension(m,m,n+1) :: nt !n_1 = n_0, ..., n_201 = n_200
double precision, intent(inout), dimension(m,n+1) :: rt !same as n, r_1 = r_0 etc.
double precision, intent(inout), dimension(m,n) :: ahat
double precision, intent(inout), dimension(m,m,n) :: vvt
double precision, intent(in), dimension(m,m,d+1) ::  pinf
double precision, intent(in), dimension(m,m,d+1) ::  pstar
double precision, intent(in),dimension(m,p,d) ::  kinfuni
double precision, intent(in),dimension(m,p,d) ::  kstaruni
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
double precision, dimension(m,m) ::  lt
double precision, dimension(m,m) :: nrec
double precision, dimension(m,m) :: nrec1
double precision, dimension(m,m) :: nrec2
double precision, dimension(m) :: rrec
double precision, dimension(m) :: rrec1
double precision, dimension(m) :: rhelp
double precision, dimension(m,m) ::  im
double precision, dimension(m,m) ::  mm
double precision, dimension(m,m) ::  mm2
double precision, intent(in) :: eps
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
double precision, dimension(p,p,(n-1)*tvh+1) :: ldl 
double precision, dimension(p,(n-1)*tvh+1) :: diag
integer, dimension(p) :: vp,apu
integer, dimension((n-1)*timevar(4)+1) :: hdiagtest
integer, dimension(n) :: ydimt

do t = 1, n
   ydimt(t) = sum(ymiss(1:p,t))
end do

do i=1,p
   vp(i)=i
end do

hdiagtest=0


do t=1, (n-1)*tvh+1
   ldl(1:p,1:p,t) = ht(1:p,1:p,(t-1)*timevar(4)+1)
   do i = 1, p
      diag(i,t) = ldl(i,i,t)
   end do
end do


hdiagtest=0

if(p>1) then
   do t= 1, ((n-1)*timevar(4)+1)
      test: do i = 1, ydimt(t)
         do k = i+1, ydimt(t)
            if(abs(ht(k,i,t)) > eps) then
               hdiagtest(t)=1
               exit test
            end if
         end do
      end do test
   end do 
   if(sum(hdiagtest)>0) then
      do t=1, ((n-1)*timevar(4)+1)
         call dpotrf('l',p,ldl(1:p,1:p,t),p,info)
         do i = 1, p
            diag(i,t) = ldl(i,i,t)
            ldl(1:p,i,t) = ldl(1:p,i,t)/diag(i,t)
         end do
         diag(1:p,t) = diag(1:p,t)**2
      end do
      do t=1,(n-1)*tvhz+1
         call dtrsm('l','l','n','u',p,m,1.0d0,ldl(1:p,1:p,(t-1)*timevar(4)+1),p,zt(1:p,1:m,t),p) !solve z*=inv(L) * zt
      end do
   end if
   if(yna==1) then
      if(tvhz==0) then !ht ja zt ei riipu ajasta
         do t=2,n
            zt(1:p,1:m,t) = zt(1:p,1:m,1)
         end do
      end if
      if(timevar(4)==0) then
         do t=2,n
            ldl(1:p,1:p,t) = ldl(1:p,1:p,1)
         end do
      end if
      do t = 1, n
         if ((ydimt(t) .NE. p) .AND.( ydimt(t) .NE. 0)) then
            apu=ymiss(1:p,t)*vp
            ldl(1:ydimt(t), 1:ydimt(t),t) = ldl(pack(apu,apu>0), pack(apu,apu>0), t)
            zt(1:ydimt(t),1:m , t) = zt(pack(apu,apu>0), 1:m, t)
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
         call dger(m,m,-1.0d0,ktuni(1:m,i,t),1,zt(i,1:m,(t-1)*tvhz+1),1,lt,m) !l = I -kz 
         call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1) 
         rrec = rhelp
         call daxpy(m,vtuni(i,t)/ftuni(i,t),zt(i,1:m,(t-1)*tvhz+1),1,rrec,1)
         call dsymm('l','u',m,m,1.0d0,nrec,m,lt,m,0.0d0,mm,m) !n*l
         call dgemm('t','n',m,m,m,1.0d0,lt,m,mm,m,0.0d0,nrec,m) !n = l'nl
         call dger(m,m,(1.0d0/ftuni(i,t)),zt(i,1:m,(t-1)*tvhz+1),1,zt(i,1:m,(t-1)*tvhz+1),1,nrec,m) ! n = n+z'z/f  
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
            lt = im
            call dger(m,m,-1.0d0,ktuni(1:m,i,t),1,zt(i,1:m,(t-1)*tvhz+1),1,lt,m) !l = i -kz
            call dgemv('t',m,m,1.0d0,lt,m,rrec,1,0.0d0,rhelp,1) 
            rrec=rhelp
            call daxpy(m,vtuni(i,t)/ftuni(i,t),zt(i,1:m,(t-1)*tvhz+1),1,rrec,1)
            call dgemm('n','n',m,m,m,1.0d0,nrec,m,lt,m,0.0d0,mm,m) !n*l
            call dgemm('t','n',m,m,m,1.0d0,lt,m,mm,m,0.0d0,nrec,m) !n = l'nl
            call dger(m,m,(1.0d0)/ftuni(i,t),zt(i,1:m,(t-1)*tvhz+1),1,zt(i,1:m,(t-1)*tvhz+1),1,nrec,m) ! n = n+z'z/f      
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
         call dger(m,m,-1.0d0/finfuni(i,t),kinfuni(1:m,i,t),1,zt(i,1:m,(t-1)*tvhz+1),1,linfuni,m) !linf
         rhelp = -kstaruni(1:m,i,t)
         call daxpy(m,fstaruni(i,t)/finfuni(i,t),kinfuni(1:m,i,t),1,rhelp,1)
         l0=0.0d0
         call dger(m,m,(1.0d0/finfuni(i,t)),rhelp,1,zt(i,1:m,(t-1)*tvhz+1),1,l0,m) !l0

         call dgemv('t',m,m,1.0d0,linfuni,m,rrec1,1,0.0d0,rhelp,1) !rt1
         call dcopy(m,rhelp,1,rrec1,1)
         call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
         call daxpy(m,(vtuni(i,t)/finfuni(i,t)),zt(i,1:m,(t-1)*tvhz+1),1,rrec1,1)

         call dgemv('t',m,m,1.0d0,linfuni,m,rrec,1,0.0d0,rhelp,1) !rt0 
         rrec = rhelp      

 
         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
         call dgemm('n','n',m,m,m,1.0d0,mm,m,linfuni,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
         
         call dger(m,m,-1.0d0*fstaruni(i,t)/(finfuni(i,t)**2.0d0),zt(i,1:m,(t-1)*tvhz+1),1,zt(i,1:m,(t-1)*tvhz+1),1,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2
         
         call dsymm('l','u',m,m,1.0d0,nrec,m,l0,m,0.0d0,mm,m) !mm= nt0*l0
         
         call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2 + l0'*nt0*l0
         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
         call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = nt2 + linf'*nt1*l0
         call dgemm('t','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1'*linf
         call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1'*linf hUOm ntrans
         ! call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec1,m,0.0d0,mm,m) !mm = l0'*nt1
         ! call dgemm('n','n',m,m,m,1.0d0,mm,m,linfuni,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1*linf 
         
         call dgemm('n','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1*linf !!!!!!!!!!
         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec1,m) !nt1 = linf'*mm
         call dger(m,m,(1.0d0)/finfuni(i,t),zt(i,1:m,(t-1)*tvhz+1),1,zt(i,1:m,(t-1)*tvhz+1),1,nrec1,m) 
         !nt1 = linf'nt1'linf + z'z/finf
         call dsymm('l','u',m,m,1.0d0,nrec,m,linfuni,m,0.0d0,mm,m) !mm= nt0*linf
         call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'nt1*linf + z'z/finf
         
         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec,m) !nt0 = linf'*mm

         !         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
         !         call dgemm('n','n',m,m,m,1.0d0,mm,m,linfuni,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
         !         call dger(m,m,(-fstaruni(i,t)/(finfuni(i,t)**2.0d0)),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec2,m) 
         !         call dsymm('l','u',m,m,1.0d0,nrec,m,l0,m,0.0d0,mm,m) !mm= nt0*l0
         
         !         call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2 + l0'*nt0*l0
         !         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
         !         call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = nt2 + linf'*nt1*l0
         !       !  call dgemm('t','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1'*linf
         !       !  call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1'*linf hUOm ntrans
         !         call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec1,m,0.0d0,mm,m) !mm = l0'*nt1
         !         call dgemm('n','n',m,m,m,1.0d0,mm,m,linfuni,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1*linf 

         !         call dgemm('n','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1*linf !!!!!!!!!!
         !         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec1,m) !nt1 = linf'*mm
         !         call dger(m,m,(1.0d0)/finfuni(i,t),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec1,m) 
         !         call dsymm('l','u',m,m,1.0d0,nrec,m,linfuni,m,0.0d0,mm,m) !mm= nt0*linf
         !         call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'*nt1*linf + z'z/finf
      
         !         call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec,m) !nt0 = linf'*mm

         
      else
         lstaruni= im
         call dger(m,m,(-1.0d0)/fstaruni(i,t),kstaruni(1:m,i,t),1,zt(i,1:m,(t-1)*tvhz+1),1,lstaruni,m) !lstar = I -Kstar*Z/Fstar         
         call dgemv('t',m,m,1.0d0,lstaruni,m,rrec,1,0.0d0,rhelp,1)
         rrec = rhelp
         call daxpy(m,vtuni(i,t)/fstaruni(i,t),zt(i,1:m,(t-1)*tvhz+1),1,rrec,1) !r0 = Z'vtuni/Fstar - Lstar'r0
         call dgemv('t',m,m,1.0d0,lstaruni,m,rrec1,1,0.0d0,rhelp,1)
         rrec1=rhelp         
         
         call dgemm('t','n',m,m,m,1.0d0,lstaruni,m,nrec,m,0.0d0,mm,m) !mm =lstar'*nt0
         call dgemm('n','n',m,m,m,1.0d0,mm,m,lstaruni,m,0.0d0,nrec,m) !nt0 = lstar'*nt0*lstar
         call dger(m,m,(1.0d0)/fstaruni(i,t),zt(i,1:m,(t-1)*tvhz+1),1,zt(i,1:m,(t-1)*tvhz+1),1,nrec,m)  !nt0 = z'z/fstar+lstar'*nt0*lstar
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
            call dger(m,m,-1.0d0/finfuni(i,t),kinfuni(1:m,i,t),1,zt(i,1:m,(t-1)*tvhz+1),1,linfuni,m) !linf
            rhelp = -kstaruni(1:m,i,t)
            call daxpy(m,fstaruni(i,t)/finfuni(i,t),kinfuni(1:m,i,t),1,rhelp,1)
            l0=0.0d0
            call dger(m,m,(1.0d0/finfuni(i,t)),rhelp,1,zt(i,1:m,(t-1)*tvhz+1),1,l0,m) !l0

            call dgemv('t',m,m,1.0d0,linfuni,m,rrec1,1,0.0d0,rhelp,1) !rt1
            call dcopy(m,rhelp,1,rrec1,1)
            call dgemv('t',m,m,1.0d0,l0,m,rrec,1,1.0d0,rrec1,1)
            call daxpy(m,vtuni(i,t)/finfuni(i,t),zt(i,1:m,(t-1)*tvhz+1),1,rrec1,1)
            call dgemv('t',m,m,1.0d0,linfuni,m,rrec,1,0.0d0,rhelp,1) !rt0 
            rrec = rhelp
            
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec2,m,0.0d0,mm,m) !mm =linf'*nt2
            call dgemm('n','n',m,m,m,1.0d0,mm,m,linfuni,m,0.0d0,nrec2,m) !nt2 = linf'*nt2*linf
            
            call dger(m,m,-1.0d0*fstaruni(i,t)/(finfuni(i,t)**2.0d0),zt(i,1:m,t),1,zt(i,1:m,t),1,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2
  
            call dsymm('l','u',m,m,1.0d0,nrec,m,l0,m,0.0d0,mm,m) !mm= nt0*l0
           
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = linf'nt2'linf + z'z*fstar/finf^2 + l0'*nt0*l0
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,nrec1,m,0.0d0,mm,m) !mm = linf'*nt1
            call dgemm('n','n',m,m,m,1.0d0,mm,m,l0,m,1.0d0,nrec2,m) !nt2 = nt2 + linf'*nt1*l0
            call dgemm('t','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1'*linf
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1'*linf hUOm ntrans
           ! call dgemm('t','n',m,m,m,1.0d0,l0,m,nrec1,m,0.0d0,mm,m) !mm = l0'*nt1
           ! call dgemm('n','n',m,m,m,1.0d0,mm,m,linfuni,m,1.0d0,nrec2,m) !nt2 = nt2 + l0'*nt1*linf 

            call dgemm('n','n',m,m,m,1.0d0,nrec1,m,linfuni,m,0.0d0,mm,m) !mm = nt1*linf !!!!!!!!!!
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec1,m) !nt1 = linf'*mm
            call dger(m,m,(1.0d0)/finfuni(i,t),zt(i,1:m,(t-1)*tvhz+1),1,zt(i,1:m,(t-1)*tvhz+1),1,nrec1,m) 
            !nt1 = linf'nt1'linf + z'z/finf
            call dsymm('l','u',m,m,1.0d0,nrec,m,linfuni,m,0.0d0,mm,m) !mm= nt0*linf
            call dgemm('t','n',m,m,m,1.0d0,l0,m,mm,m,1.0d0,nrec1,m) !nt1 = l0'*nt0*linf+ linf'nt1*linf + z'z/finf
            
            call dgemm('t','n',m,m,m,1.0d0,linfuni,m,mm,m,0.0d0,nrec,m) !nt0 = linf'*mm
            
           
         else
            lstaruni= im
            call dger(m,m,(-1.0d0)/fstaruni(i,t),kstaruni(1:m,i,t),1,zt(i,1:m,(t-1)*tvhz+1),1,lstaruni,m) !lstar = I -Kstar*Z/Fstar         
            call dgemv('t',m,m,1.0d0,lstaruni,m,rrec,1,0.0d0,rhelp,1) !oli beta 1.0d0!!!!... JA miinusmerkki
            rrec = rhelp
            call daxpy(m,vtuni(i,t)/fstaruni(i,t),zt(i,1:m,(t-1)*tvhz+1),1,rrec,1) !r0 = Z'vtuni/Fstar - Lstar'r0
            call dgemv('t',m,m,1.0d0,lstaruni,m,rrec1,1,0.0d0,rhelp,1)
            rrec1=rhelp           
            
            call dgemm('t','n',m,m,m,1.0d0,lstaruni,m,nrec,m,0.0d0,mm,m) !mm =lstar'*nt0
            call dgemm('n','n',m,m,m,1.0d0,mm,m,lstaruni,m,0.0d0,nrec,m) !nt0 = lstar'*nt0*lstar
            call dger(m,m,(1.0d0)/fstaruni(i,t),zt(i,1:m,(t-1)*tvhz+1),1,zt(i,1:m,(t-1)*tvhz+1),1,nrec,m)  !nt0 = z'z/fstar+lstar'*nt0*lstar
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

end subroutine
