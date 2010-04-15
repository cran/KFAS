subroutine kf(yt, ydimt, timevar, zt, tt, rt, ht, qt, a1, p1, at, pt, vtuni,&
     ftuni, ktuni, pinf, pstar, finfuni, fstaruni, kinfuni, kstaruni, d, j, p, m, r, n,&
     lik, optcal, info, vt, ft, kt, lt, finf, fstar, kinf, kstar, linf, lstar, eps)

implicit none

integer, intent(in) ::  p, m, r, n
integer, intent(inout) :: d, j
integer, intent(inout), dimension(4) :: info
integer ::  t, i, k
integer, intent(in), dimension(n) :: ydimt
integer, intent(in), dimension(5) :: timevar
double precision, intent(inout), dimension(p,n) :: yt
double precision, intent(inout), dimension(p,m,(n-1)*timevar(5)+1) :: zt  
double precision, intent(in), dimension(m,m,(n-1)*timevar(1)+1) :: tt 
double precision, intent(in), dimension(m,r,(n-1)*timevar(2)+1) :: rt 
double precision, intent(inout), dimension(p,p,(n-1)*timevar(4)+1) :: ht 
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
double precision, intent(inout) :: lik
double precision, dimension(p,n) :: yy
double precision, dimension(p,m,(n-1)*timevar(5)+1) :: zz 
double precision, dimension(p,p,(n-1)*timevar(4)+1) :: hh
integer, intent(in), dimension(4) :: optcal 
double precision, dimension(m) :: arec
double precision, dimension(m,m) :: prec
double precision, dimension(m,m) ::  psrec
double precision, dimension(m,m) ::  pirec
double precision, dimension(m,r) :: mr
double precision, dimension(p,m) :: pm
double precision, dimension(m,m) ::  im
double precision, dimension(m) ::  m1
double precision, dimension(m,m) ::  mm
double precision, dimension(m,p) ::  mp
double precision, dimension(p,p) ::  cholft
integer, dimension(n) :: hdiagtest
double precision, intent(in) :: eps
double precision, dimension(p) :: diag
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

if(sum(optcal)>0) then
   yy=yt
   zz=zt
   hh=ht
end if

hdiagtest = 0

if(p>1) then !testing diagonality
   do t= 1, (n-1)*timevar(4)+1
      test: do i = 1, ydimt(t)
         do k = i+1, ydimt(t)
            if(abs(ht(k,i,t)) > eps) then
               hdiagtest(t)=1
               exit test
            end if
         end do
      end do test
   end do
   if(sum(hdiagtest)/=0) then 
      do t = 1, (n-1)*timevar(4)+1
         if(hdiagtest(t)==1 .AND. ydimt(t)>0) then           
            call dpotrf('l',ydimt(t),ht(1:ydimt(t),1:ydimt(t),t),ydimt(t),info(1))
            if(info(1) /=0) then
               info(1)=1
               return
            end if
            do i = 1, ydimt(t)
               diag(i)=ht(i,i,t)
               ht(1:ydimt(t),i,t) =  ht(1:ydimt(t),i,t)/diag(i)
            end do
            do i = 1, ydimt(t)
               ht(i,i,t) = diag(i)**2 
            end do
         end if
      end do
      do t = 1,(n-1)*max(timevar(5),timevar(4))+1
         if(ydimt(t)>0 .AND. hdiagtest(t)/=0) then
            call dtrsm('l','l','n','u',ydimt(t),m,1.0d0,ht(1:ydimt(t),1:ydimt(t),(t-1)*timevar(4)+1)&
                 ,ydimt(t),zt(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t)) !solve z*=inv(L) * zt
         end if
      end do
      do t= 1, n
         if(ydimt(t)>0) then           
            call dtrsv('l','n','u',ydimt(t),ht(1:ydimt(t),1:ydimt(t),(t-1)*timevar(4)+1),ydimt(t),yt(1:ydimt(t),t),1) !solve y*=inv(L) * yt, 
         end if
      end do
   end if
end if

lik = 0.0d0

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
         vtuni(j,d) = yt(j,d) - ddot(m,zt(j,1:m,(d-1)*timevar(5)+1),1,arec,1) !arec
         call dsymv('u',m,1.0d0,psrec,m,zt(j,1:m,(d-1)*timevar(5)+1),1,0.0d0,kstaruni(1:m,j,d),1) ! kstar_t,i = pstar_t,i*t(z_t,i)
         call dsymv('u',m,1.0d0,pirec,m,zt(j,1:m,(d-1)*timevar(5)+1),1,0.0d0,kinfuni(1:m,j,d),1) ! kinf_t,i = pinf_t,i*t(z_t,i)
         if (abs(finfuni(j,d)) > eps) then
            call daxpy(m,vtuni(j,d)/finfuni(j,d),kinfuni(1:m,j,d),1,arec,1) !a_rec = a_rec + kinf(:,i,t)*vt(:,t)/finf(j,d)
            call dsyr('u',m,fstaruni(j,d)/(finfuni(j,d)**2),kinfuni(1:m,j,d),1,psrec,m) !psrec = psrec +  kinf*kinf'*fstar/finf^2
            call dsyr2('u',m,-1.0d0/finfuni(j,d),kstaruni(1:m,j,d),1,kinfuni(1:m,j,d),1,psrec,m) !psrec = psrec -(kstar*kinf'+kinf*kstar')/finf
            
            !call dger(m,m,fstaruni(j,d)/(finfuni(j,d)**2),kinfuni(1:m,j,d),1,kinfuni(1:m,j,d),1,psrec,m) !psrec = psrec +  kinf*kinf'*fstar/finf^2
            !call dger(m,m,-1.0d0/finfuni(j,d),kstaruni(1:m,j,d),1,kinfuni(1:m,j,d),1,psrec,m)
            !call dger(m,m,-1.0d0/finfuni(j,d),kinfuni(1:m,j,d),1,kstaruni(1:m,j,d),1,psrec,m)
            
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
            call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(d-1)*timevar(3)+1),r,rt(1:m,1:r,(d-1)*timevar(2)+1),m,0.0d0,mr,m)
            call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(1:m,1:r,(d-1)*timevar(2)+1),m,1.0d0,pstar(1:m,1:m,d+1),m)   
         else
            call dger(m,m,qt(1,1,(d-1)*timevar(3)+1),rt(1:m,1,(d-1)*timevar(2)+1),1,rt(1:m,1,(d-1)*timevar(2)+1),1,&
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
      vtuni(i,d) = yt(i,d) - ddot(m,zt(i,1:m,(d-1)*timevar(5)+1),1,arec,1) !vtuni         
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
         call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(d-1)*timevar(3)+1),r,rt(1:m,1:r,(d-1)*timevar(2)+1),m,0.0d0,mr,m)
         call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(1:m,1:r,(d-1)*timevar(2)+1),m,1.0d0,pt(1:m,1:m,d+1),m)   
      else
         call dger(m,m,qt(1,1,(d-1)*timevar(3)+1),rt(1:m,1,(d-1)*timevar(2)+1),1,rt(1:m,1,(d-1)*timevar(2)+1),1,&
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
      vtuni(i,t) = yt(i,t) - ddot(m,zt(i,1:m,(t-1)*timevar(5)+1),1,arec,1) !univariate vt           
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
   !pt(:,:,t+1) = matmul(matmul(tt,p_rec),transpose(tt)) + matmul(matmul(rt,qt),transpose(rt))
   call dsymm('r','u',m,m,1.0d0,prec,m,tt(1:m,1:m,(t-1)*timevar(1)+1),m,0.0d0,mm,m)
   call dgemm('n','t',m,m,m,1.0d0,mm,m,tt(1:m,1:m,(t-1)*timevar(1)+1),m,0.0d0,pt(1:m,1:m,t+1),m)
   if(m/=r) then
      if(r>1) then
         call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,rt(1:m,1:r,(t-1)*timevar(2)+1),m,0.0d0,mr,m)
         call dgemm('n','t',m,m,r,1.0d0,mr,m,rt(1:m,1:r,(t-1)*timevar(2)+1),m,1.0d0,pt(1:m,1:m,t+1),m)   
      else
         call dger(m,m,qt(1,1,(t-1)*timevar(3)+1),rt(1:m,1,(t-1)*timevar(2)+1),1,rt(1:m,1,(t-1)*timevar(2)+1),& 
              1,pt(1:m,1:m,t+1),m) 
      end if
   else
      pt(1:m,1:m,t+1) = pt(1:m,1:m,t+1) + qt(1:r,1:r,(t-1)*timevar(3)+1)
   end if
   call dcopy(m,at(1:m,t+1),1,arec,1) ! a_rec =at(:,t+1)
   prec = pt(1:m,1:m,t+1)
end do


!real ft, vt, kt, lt:

if(sum(optcal)>0) then

if(optcal(1)==1) then
   do t=1, n
      if(ydimt(t)>0) then
         if(ydimt(t)==1) then
            vt(1,t)=vtuni(1,t)
         else
            call dcopy(ydimt(t),yy(1:ydimt(t),t),1,vt(1:ydimt(t),t),1)
            call dgemv('n',ydimt(t),m,-1.0d0,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t),at(1:m,t),1,1.0d0,vt(1:ydimt(t),t),1)
         end if
      end if
   end do
end if

if(optcal(2)==1) then
   do t= 1, d
      if(ydimt(t)>0) then
         if(p==1) then
            fstar(1,1,t) = fstaruni(1,t)
            finf(1,1,t) = finfuni(1,t)
         else
            call dgemm('n','t',m,ydimt(t),m,1.0d0,pinf(1:m,1:m,t),m,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t),&
                 0.0d0,mp(1:m,1:ydimt(t)),m) !mp = pinf*z'
            call dgemm('n','n',ydimt(t),ydimt(t),m,1.0d0,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),p,&
                 mp(1:m,1:ydimt(t)),m,0.0d0,finf(1:ydimt(t),1:ydimt(t),t),ydimt(t)) !finf = z*mp
            
            call dgemm('n','t',m,ydimt(t),m,1.0d0,pstar(1:m,1:m,t),m,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),&
                 ydimt(t),0.0d0,mp(1:m,1:ydimt(t)),m) !mp = pstar*z'
            call dgemm('n','n',ydimt(t),ydimt(t),m,1.0d0,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t),&
                 mp(1:m,1:ydimt(t)),m,0.0d0,fstar(1:ydimt(t),1:ydimt(t),t),ydimt(t)) !fstar = z*mp 
            fstar(1:ydimt(t),1:ydimt(t),t) = fstar(1:ydimt(t),1:ydimt(t),t) +  hh(1:ydimt(t),1:ydimt(t),(t-1)*timevar(4)+1)
         end if
      end if
   end do
   do t=d+1, n
      if(ydimt(t)>0) then
         if(p==1) then
            ft(1,1,t) = ftuni(1,t)
         else
            ft(1:ydimt(t),1:ydimt(t),t) =  hh(1:ydimt(t),1:ydimt(t),(t-1)*timevar(4)+1)
            call dsymm('r','u',ydimt(t),m,1.0d0,pt(1:m,1:m,t),m,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t),&
                 0.0d0,pm(1:ydimt(t),1:m),ydimt(t))
            call dgemm('n','t',ydimt(t),ydimt(t),m,1.0d0,pm(1:ydimt(t),1:m),ydimt(t),&
                 zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t),1.0d0,ft(1:ydimt(t),1:ydimt(t),t),ydimt(t)) 
         end if
      end if
   end do
end if


if(optcal(3)==1 .AND. optcal(2)==1) then   
   mulkd:  do t=1,d
      if(ydimt(t)>0) then
         if(maxval(abs(finf(1:ydimt(t),1:ydimt(t),t)))<eps) then !finf=0
            call dgemm('n','t',m,ydimt(t),m,1.0d0,pstar(1:m,1:m,t),m,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),p,0.0d0,&
                 mp(1:m,1:ydimt(t)),m) !mp = pstar*z'
            call dgemm('n','n',m,ydimt(t),m,1.0d0,tt(1:m,1:m,(t-1)*timevar(1)+1),m,mp(1:m,1:ydimt(t)),m,0.0d0,&
                 kstar(1:m,1:ydimt(t),t),m) !kstar = t*mp
            pm(1:ydimt(t),1:m) = transpose(kstar(1:m,1:ydimt(t),t))
            cholft(1:ydimt(t),1:ydimt(t)) = fstar(1:ydimt(t),1:ydimt(t),t)
            call dposv('l',ydimt(t),m,cholft(1:ydimt(t),1:ydimt(t)),ydimt(t),pm(1:ydimt(t),1:m),ydimt(t),info(2))
            if(info(2) /=0) then
               info(2)=1
               exit mulkd
            end if
            kstar(1:m,1:ydimt(t),t) = transpose(pm(1:ydimt(t),1:m))
         else
            pm(1:ydimt(t),1:m) = zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1)
            cholft(1:ydimt(t),1:ydimt(t)) = finf(1:ydimt(t),1:ydimt(t),t)
            call dposv('l',ydimt(t),m,cholft(1:ydimt(t),1:ydimt(t)),ydimt(t),pm(1:ydimt(t),1:m),ydimt(t),info(3))
            if(info(3)/=0) then !if finf not pos.def
               info(3)=1
               exit mulkd
            end if
            kinf(1:m,1:ydimt(t),t) = transpose(pm(1:ydimt(t),1:m)) !kinf = z'*inv(finf)
            
            call dgemm('n','n',m,ydimt(t),m,1.0d0,pinf(1:m,1:m,t),m,kinf(1:m,1:ydimt(t),t),m,0.0d0,&
                 mp(1:m,1:ydimt(t)),m) !mp = pinf*kinf*z'*inv(finf)
            call dgemm('n','n',m,ydimt(t),m,1.0d0,tt(1:m,1:m,(t-1)*timevar(1)+1),m,mp(1:m,1:ydimt(t)),&
                 m,0.0d0,kinf(1:m,1:ydimt(t),t),m) !kinf = t*pinf*z'*inv(finf)
               
            call dgemm('n','t',m,ydimt(t),m,1.0d0,pstar(1:m,1:m,t),m,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t),0.0d0,&
                    mp(1:m,1:ydimt(t)),m) !mp = pstar*z'    
            call dgemm('n','n',m,ydimt(t),m,1.0d0,tt(1:m,1:m,(t-1)*timevar(1)+1),m,mp(1:m,1:ydimt(t)),m,&
                 0.0d0,kstar(1:m,1:ydimt(t),t),m) !kstar = t*mp
            
            call dgemm('n','n',m,ydimt(t),ydimt(t),1.0d0,kinf(1:m,1:ydimt(t),t),m,&
                 fstar(1:ydimt(t),1:ydimt(t),t),ydimt(t),0.0d0,mp(1:m,1:ydimt(t)),m) !mp =kinf*fstar
            kstar(1:m,1:ydimt(t),t) = kstar(1:m,1:ydimt(t),t) - mp(1:m,1:ydimt(t)) ! note the -, wrong in formula of Dk2003, corrected form in appendix!
            pm(1:ydimt(t),1:m) = transpose(kstar(1:m,1:ydimt(t),t))
            call dpotrs('l',ydimt(t),m,cholft(1:ydimt(t),1:ydimt(t)),ydimt(t),pm(1:ydimt(t),1:m),ydimt(t),info(3)) !use cholesky from earlier point               
            if(info(3) /=0) then
               info(3)=1
                exit mulkd
            end if
            kstar(1:m,1:ydimt(t),t) = transpose(pm(1:ydimt(t),1:m))     
         end if
      end if
   end do mulkd
   mulk: do t=d+1, n
      if(ydimt(t)>0) then
         if(p==1) then
            call dgemv('n',m,m,1.0d0,tt(1:m,1:m,(t-1)*timevar(1)+1),m,ktuni(1:m,1,t),1,0.0d0,kt(1:m,1,t),1) !kt=t*ktuni
         else
            call dsymm('r','u',m,m,1.0d0,pt(1:m,1:m,t),m,tt(1:m,1:m,(t-1)*timevar(1)+1),m,0.0d0,mm,m) !TP
            call dgemm('n','t',m,ydimt(t),m,1.0d0,mm,m,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t),&
                 1.0d0,kt(1:m,1:ydimt(t),t),m) !TPZ'
            cholft(1:ydimt(t),1:ydimt(t)) = transpose(ft(1:ydimt(t),1:ydimt(t),t))
            pm(1:ydimt(t),1:m) = transpose(kt(1:m,1:ydimt(t),t))
            call dposv('l',ydimt(t),m,cholft(1:ydimt(t),1:ydimt(t)),ydimt(t),pm(1:ydimt(t),1:m),ydimt(t),info(4))
            if(info(4) /=0) then
               info(4)=1
               exit mulk
            end if
            kt(1:m,1:ydimt(t),t) = transpose(pm(1:ydimt(t),1:m))        
         end if
      end if
   end do mulk
   
end if

if(optcal(4)==1 .AND. optcal(3)==1  .AND. sum(info)==0) then
   do t=1, d
      if(ydimt(t)>0) then
         if(maxval(abs(finf(1:ydimt(t),1:ydimt(t),t)))<eps) then !finf=0
            lstar(1:m,1:m,t) = tt(1:m,1:m,(t-1)*timevar(1)+1)
            call dgemm('n','n',m,m,ydimt(t),-1.0d0,kstar(1:m,1:ydimt(t),t),m,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t),&
                 1.0d0,lstar(1:m,1:m,t),m) !lstar
         else            
            linf(1:m,1:m,t) = tt(1:m,1:m,(t-1)*timevar(1)+1)
            call dgemm('n','n',m,m,ydimt(t),-1.0d0,kinf(1:m,1:ydimt(t),t),m,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t),&
                 1.0d0,linf(1:m,1:m,t),m) !linf 
            call dgemm('n','n',m,m,ydimt(t),-1.0d0,kstar(1:m,1:ydimt(t),t),m,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1),ydimt(t),&
                 0.0d0,lstar(1:m,1:m,t),m) !lstar
         end if
      end if
   end do
   do t = d+1, n
      lt(1:m,1:m,t) = tt(1:m,1:m,(t-1)*timevar(1)+1)
      if(ydimt(t)>0) then
         call dgemm('n','n',m,m,ydimt(t),-1.0d0,kt(1:m,1:ydimt(t),t),m,zz(1:ydimt(t),1:m,(t-1)*timevar(5)+1)&
              ,ydimt(t),1.0d0,lt(1:m,1:m,t),m) !lt = t - kz
      end if
   end do
end if

end if

end subroutine
 
