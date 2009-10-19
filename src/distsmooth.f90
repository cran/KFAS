subroutine distsmooth(timevar, ht,rtv,qt,ft,kt,vt,nt,rt,fstar,finf,&
kinf,kstar,nt0,rt0,d,p,m,r,n,eps,epshat,epshatvar,etahat,etahatvar,info)

implicit none
integer, intent(in) :: d, p, m, r, n
integer :: t
integer, intent(in), dimension(3) :: timevar
double precision, intent(in) :: eps
double precision, intent(in), dimension(p,p,(n-1)*timevar(1)+1) :: ht
double precision, intent(in), dimension(m,r,(n-1)*timevar(2)+1) :: rtv
double precision, intent(in), dimension(r,r,(n-1)*timevar(3)+1) :: qt

double precision, intent(inout), dimension(p,n) ::  vt
double precision, intent(inout), dimension(p,p,n) :: ft
double precision, intent(in), dimension(m,p,n) :: kt

double precision, intent(in),dimension(m,p,d) ::  kinf
double precision, intent(in),dimension(m,p,d) ::  kstar
double precision, intent(inout),dimension(p,p,d) ::  fstar
double precision, intent(in),dimension(p,p,d) ::  finf

double precision, intent(in), dimension(m,m,n+1) :: nt 
double precision, intent(in), dimension(m,n+1) :: rt
double precision, intent(in), dimension(m,d+1) :: rt0
double precision, intent(in), dimension(m,m,d+1) :: nt0

double precision, intent(inout), dimension(p,n) :: epshat
double precision, intent(inout), dimension(p,p,n) :: epshatvar
double precision, intent(inout), dimension(r,n) :: etahat
double precision, intent(inout), dimension(r,r,n) :: etahatvar
double precision, dimension(m) :: help
double precision, dimension(m,r) :: mr
double precision, dimension(m,r) :: mr2
double precision, dimension(m,p) :: mp
double precision, dimension(m,p) :: mp2
double precision, dimension(p,p) :: pp
integer, intent(inout) ::  info

epshatvar = 0.0d0
do t= 1, p
   epshatvar(t,t,1:n) = 1.0d0
end do

do t = 1, d
   if(maxval(abs(finf(1:p,1:p,t))) < eps) then
      call dposv('l',p,1,fstar(1:p,1:p,t),p,vt(1:p,t),p,info)
      if(info /=0) then
         info=1
         return
      end if
      call dgemv('t',m,p,-1.0d0,kstar(1:m,1:p,t),m,rt0(1:m,t+1),1,1.0d0,vt(1:p,t),1)
      call dsymv('l',p,1.0d0,ht(1:p,1:p,(t-1)*timevar(1)+1),p,vt(1:p,t),1,0.0d0,epshat(1:p,t),1)

      call dpotrs('l',p,p,fstar(1:p,1:p,t),p,epshatvar(1:p,1:p,t),p,info)
      call dgemm('n','n',m,p,m,1.0d0,nt0(1:m,1:m,t+1),m,kstar(1:m,1:p,t),m,0.0d0,mp,m)
      call dgemm('t','n',p,p,m,1.0d0,kstar(1:m,1:p,t),m,mp(1:m,1:p),m,1.0d0,epshatvar(1:p,1:p,t),p)
      call dsymm('l','u',p,p,1.0d0,ht(1:p,1:p,(t-1)*timevar(1)+1),p,epshatvar(1:p,1:p,t),p,0.0d0,pp,p)
      call dsymm('r','u',p,p,1.0d0,ht(1:p,1:p,(t-1)*timevar(1)+1),p,pp,p,0.0d0,epshatvar(1:p,1:p,t),p)
      epshatvar(1:p,1:p,t) = ht(1:p,1:p,(t-1)*timevar(1)+1) - epshatvar(1:p,1:p,t)

      call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,rt0(1:m,t+1),1,0.0d0,help,1)
      call dsymv('l',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,help,1,0.0d0,etahat(1:r,t),1)
      etahatvar(1:r,1:r,t) = qt(1:r,1:r,(t-1)*timevar(3)+1)
      call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,0.0d0,mr,m)
      call dgemm('n','n',m,r,m,1.0d0,nt0(1:m,1:m,t+1),m,mr,m,0.0d0,mr2,m)
      call dgemm('t','n',r,r,m,-1.0d0,mr,m,mr2,m,1.0d0,etahatvar(1:r,1:r,t),r)
   else
      call dgemv('t',m,p,1.0d0,kinf(1:m,1:p,t),m,rt0(1:m,t+1),1,0.0d0,help,1)
      call dsymv('l',p,-1.0d0,ht(1:p,1:p,(t-1)*timevar(1)+1),p,help,1,0.0d0,epshat(1:p,t),1)      
      epshatvar(1:p,1:p,t) = ht(1:p,1:p,(t-1)*timevar(1)+1)
      call dsymm('r','u',m,p,1.0d0,ht(1:p,1:p,(t-1)*timevar(1)+1),p,kinf(1:m,1:p,t),m,0.0d0,mp,m)
      call dgemm('n','n',m,p,m,1.0d0,nt0(1:m,1:m,t+1),m,mp,m,0.0d0,mp2,m)
      call dgemm('t','n',p,p,m,-1.0d0,mp,m,mp2,m,1.0d0,epshatvar(1:p,1:p,t),p)

      call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,rt0(1:m,t+1),1,0.0d0,help,1)
      call dsymv('l',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,help,1,0.0d0,etahat(1:r,t),1)
      etahatvar(1:r,1:r,t) = qt(1:r,1:r,(t-1)*timevar(3)+1)
      call dsymm('r','u',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,0.0d0,mr,m)
      call dgemm('n','n',m,r,m,1.0d0,nt0(1:m,1:m,t+1),m,mr,m,0.0d0,mr2,m)
      call dgemm('t','n',r,r,m,-1.0d0,mr,m,mr2,m,1.0d0,etahatvar(1:r,1:r,t),r)
   end if
end do


do t= d+1,n
   call dposv('l',p,1,ft(1:p,1:p,t),p,vt(1:p,t),p,info)
   if(info /=0) then
      info=2
      return
   end if
   call dgemv('t',m,p,-1.0d0,kt(1:m,1:p,t),m,rt(1:m,t+1),1,1.0d0,vt(1:p,t),1)
   call dsymv('l',p,1.0d0,ht(1:p,1:p,(t-1)*timevar(1)+1),p,vt(1:p,t),1,0.0d0,epshat(1:p,t),1)
  
   call dpotrs('l',p,p,ft(1:p,1:p,t),p,epshatvar(1:p,1:p,t),p,info)
   call dgemm('n','n',m,p,m,1.0d0,nt(1:m,1:m,t+1),m,kt(1:m,1:p,t),m,0.0d0,mp,m)
   call dgemm('t','n',p,p,m,1.0d0,kt(1:m,1:p,t),m,mp(1:m,1:p),m,1.0d0,epshatvar(1:p,1:p,t),p)
   call dsymm('l','l',p,p,1.0d0,ht(1:p,1:p,(t-1)*timevar(1)+1),p,epshatvar(1:p,1:p,t),p,0.0d0,pp,p)
   call dsymm('r','l',p,p,1.0d0,ht(1:p,1:p,(t-1)*timevar(1)+1),p,pp,p,0.0d0,epshatvar(1:p,1:p,t),p)
   epshatvar(1:p,1:p,t) = ht(1:p,1:p,(t-1)*timevar(1)+1) - epshatvar(1:p,1:p,t)
   
   call dgemv('t',m,r,1.0d0,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,rt(1:m,t+1),1,0.0d0,help,1)
   call dsymv('l',r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,help,1,0.0d0,etahat(1:r,t),1)
   etahatvar(1:r,1:r,t) = qt(1:r,1:r,(t-1)*timevar(3)+1)
   call dsymm('r','l',m,r,1.0d0,qt(1:r,1:r,(t-1)*timevar(3)+1),r,rtv(1:m,1:r,(t-1)*timevar(2)+1),m,0.0d0,mr,m)
   call dgemm('n','n',m,r,m,1.0d0,nt(1:m,1:m,t+1),m,mr,m,0.0d0,mr2,m)
   call dgemm('t','n',r,r,m,-1.0d0,mr,m,mr2,m,1.0d0,etahatvar(1:r,1:r,t),r)
end do

end subroutine distsmooth
