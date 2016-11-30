! -*- Fortran -*- 77 file relaxing.fpp
!**************************************************************************
! File: relaxing.fpp
! started Sept 3, 1994 by Ue-Li Pen, upen@astro.princeton.edu
! Purpose: to implement a single step of the relaxing TVD scheme in
! curvilinear coordinates.
!
!
! Feb 6, 1995: everything seems to work well, but we will now add the
! exact energy conservation expression.  This will require three extra
! arrays to implement to second order accuracy.
!     
! Feb 11: energy conservation seems to work.  We underestimate \sigma
!         and my current suspicion is that it stems from the TVD limiter.
!     
!     
!**************************************************************************
! Note that ng3 must be 2**n >= 16 for the relaxation to work.
!     
!     
subroutine relaxing(u,def,defp,cmax,dt,nfluidcmp)
!
! we need to pass nfluidcmp as an argument, since I want to be able
! to simulate several phases if necessary.  E.g., one anisotropic dark
! matter plus one ideal gas.
! 'c' is the freezing soundspeed.  Since we use units where dx=1,
!  It must satisfy dt < 1/c.
! To calculate c, one needs to construct the whole metric, which
! is a sacrifice of computational efficiency over memory.  It should
! be a small effect compared to the time integrator.
!
implicit none
#include 'relaxgl.f90'
integer nfluidcmp
!
! both def and defp are defined at the half time step between the
! current u and the next one that this routines computes.
!
real u(nfluidcmp,ng1,ng2,ng3), def(ng1,ng2,ng3), defp(ng1,ng2,ng3)
!mydist u(*,*,block,*), def(*,block,*)
!mydist  defp(*,block,*)
real cmax,dt
!
! Notes:
!    defp must be defined for the middle of the time step, so that the
! runge kutta has grid velocities known to second order.      
!
! locals:
!
! the parameter nstore is the number of tiers of metric and various
! other temporary stores required to achieve the second order
! runge-kutta time integrator.
!
!
! arrays:
! we assume everything is statically allocated.
!
real baj(7,ng1,ng2,nstore), vt(ng1,ng2,nstore,3,maxfluidcmp),&
     res(ng1,ng2,maxfluidcmp), u1(maxfluidcmp,ng1,ng2,nstore),&
     vtf(ng1,ng2,nstore,3,maxfluidcmp)
! we need to retain a copy of the top layer:
real vtfstore(ng1,ng2,4,3,maxfluidcmp),c(ng1,ng2,ncstore),u1store(maxfluidcmp,ng1,ng2,4)
!mydist baj(*,*,block,*),vt(*,block,*,*,*),res(*,block,*)
!mydist u1(*,*,block,*),vtf(*,block,*,*,*)
!mydist vtfstore(*,block,*,*,*)
!mydist c(*,block,*),u1store(*,*,block,*)
integer indx(ng3), indxu(ng3), indxc(ng3)
!
! array indx has image range (1:nstore), and maps the desired index
! onto the shorter temporary arrays.
!
! procedure scope:
integer idxptr, idxcptr
!
! short range locals:
integer k,kg,kh,kf,j,i,nf,idim
!

!----------------------------------------------------------------------
! begin executables: initialization
!      

do k=1,ng3
   indxu(k)=k
enddo
!
! begin main body
!
idxptr=0

! build the whole c-index pointers.  One could have done the same with
! the u and vt arrays.      
idxcptr=0
! we will need the first six again at the end, and 7 running
do k=1,ncstore
   indxc(k)=k
enddo
idxcptr=6
do k=ncstore+1,ng3
   idxcptr=mod(idxcptr-6,ncstore-6)+7
   indxc(k)=idxcptr
enddo

! bootstrap the whole process:
! we need the metric over a long range:
do kg=1,4
   idxptr=mod(idxptr,nstore)+1
   indx(kg)=idxptr
   call calcbaj(baj,def,indx,kg)
   call cflkeuler(c,u,baj,defp,indxc,indx,kg)
   call flux(vt,defp,u,baj,indx,indxu,nfluidcmp,ng3,kg)
enddo

!fpp$ skip
do k=3,6
! dpv calculates the upwind TVD limited flux 
! the domain of dependence of dpv is 2 on each side.
   kg=k+2
   kh=k
   idxptr=mod(idxptr,nstore)+1
   indx(kg)=idxptr
   call calcbaj(baj,def,indx,kg)
   call cflkeuler(c,u,baj,defp,indxc,indx,kg)
   call flux(vt,defp,u,baj,indx,indxu,nfluidcmp,ng3,kg)
   call dpv(res,u,vt,c,indx,indxu,indxc,nfluidcmp,ng3,k)
! Runge-kutta half step
   call stepu(u1,u,res,dt/2,indx,nstore,nfluidcmp,kh)
   call flux(vtf,defp,u1,baj,indx,indx,nfluidcmp,nstore,kh)
enddo
! we will need these 4 vtf layers again at the very end:

!$omp parallel default(private) & 
!$omp firstprivate(nfluidcmp) &
!$omp shared(vtfstore,vtf,indx,u1store,u1)
!fpp$ skip
!c*$* assert do(serial)
do nf=1,nfluidcmp
!fpp$ skip
   do idim=1,3
!fpp$ skip

      do k=1,4
!$dir prefer_parallel_ext
!*$* assert do (concurrent)
!$omp do
         do j=1,ng2
            do i=1,ng1
               vtfstore(i,j,k,idim,nf)=vtf(i,j,indx(k+2),idim,nf)
            enddo
         enddo
!$omp enddo nowait
      enddo
   enddo
enddo
!fpp$ skip
!c*$* assert do(serial)
!cdir novector
do nf=1,nfluidcmp
   do k=1,4
!$dir prefer_parallel_ext               
!*$* assert do (concurrent)
!$omp do
      do j=1,ng2
!cdir vector
         do i=1,ng1
            u1store(nf,i,j,k)=u1(nf,i,j,indx(k+2))
         enddo
      enddo
!$omp enddo nowait
   enddo
enddo
!$omp end parallel

! main loop down the layers.  k=5 will be the first finished layer.
do k=7,ng3+2
! the main iteration variable kg tracks the level at which we update
! the gravity
   kh=mod(k+ng3-1,ng3)+1
   kg=mod(kh+1,ng3)+1
   idxptr=mod(idxptr,nstore)+1
   indx(kg)=idxptr
! kh tracks the half step runge kutta update
   kf=mod(kh-3+ng3,ng3)+1
! kf tracks the full step runge kutta index.  At the end of the loop,
! layer kf is completely updated.
   call calcbaj(baj,def,indx,kg)
!
! flux() uses only the metric and u() at the current level, and returns
! the 3 x nfluid flux matrix at the current level.
   call cflkeuler(c,u,baj,defp,indxc,indx,kg)
   call flux(vt,defp,u,baj,indx,indxu,nfluidcmp,ng3,kg)
!
! dpv() uses u and vt at two levels each way.  It solves the system of
! constant coefficient equations
!
!  \dot{u}   + c (v^1,x+v^2,y) = 0
!  \dot{v^1} + c  u,x          = 0   ; -> split (F1-v^1)/\epsilon
!  \dot{v^2} + c        u,y    = 0   ;          (F2-v^2)/\epsilon
!
   call dpv(res,u,vt,c,indx,indxu,indxc,nfluidcmp,ng3,kh)
! Runge-Kutta half step
!
! update first order half step estimate:  u^(1/2)=u-dpv*dt/2
   call stepu(u1,u,res,dt/2,indx,nstore,nfluidcmp,kh)
! now redo to second order.
   call flux(vtf,defp,u1,baj,indx,indx,nfluidcmp,nstore,kh)
!         write(*,'(5f15.9)')(vtf(1,1,indx(kf),3,j),j=1,5)
   call dpv(res,u1,vtf,c,indx,indx,indxc,nfluidcmp,nstore,kf)
! FULL step         
! u^+ = u-dpv(u1)*dt
!         call stepu(u,u,res,dt,indxu,ng3,nfluidcmp,kf)
   call stepuu(u,res,dt,nfluidcmp,kf)
enddo

! the last four layers need to use the buffer zone:
!fpp$ skip
do kf=1,4
   kg=mod(kf+3,ng3)+1
   idxptr=mod(idxptr,nstore)+1
   indx(kg)=idxptr
!$omp parallel default(private) firstprivate(kf,nfluidcmp) & 
!$omp shared(indx,vtf,vtfstore,u1,u1store)
!fpp$ skip
!c*$* assert do prefer(serial)
   do nf=1,nfluidcmp
!fpp$ skip
      do idim=1,3

!$dir prefer_parallel_ext
!*$* assert do (concurrent)
!$omp do
         do j=1,ng2
            do i=1,ng1
             vtf(i,j,indx(kf+2),idim,nf)=vtfstore(i,j,kf,idim,nf)
            enddo
         enddo
!$omp enddo nowait
      enddo
   enddo
!fpp$ skip
!c*$* assert do(serial)
!cdir novector
   do nf=1,nfluidcmp
!$dir prefer_parallel_ext               
!*$* assert do (concurrent)
!$omp do
      do j=1,ng2
!cdir vector
         do i=1,ng1
            u1(nf,i,j,indx(kf+2))=u1store(nf,i,j,kf)
         enddo
      enddo
!$omp enddo nowait
   enddo
!$omp end parallel

   call dpv(res,u1,vtf,c,indx,indx,indxc,nfluidcmp,nstore,kf)
! FULL step         
   call stepuu(u,res,dt,nfluidcmp,kf)
enddo
      
! that was simple!  we are done.
!
return
endsubroutine relaxing


!**********************************************************************
subroutine dpv(res,u,vt,c,indxr, indxu, indxc, nfluid, ngu3, k)
implicit none
#include 'relaxgl.f90'
integer k, indxu(ng3), indxr(ng3), indxc(ng3), ngu3,nfluid
real res(ng1,ng2,nfluid), u(nfluid,ng1,ng2,ngu3),&
     vt(ng1,ng2,nstore,3,nfluid), c(ng1,ng2,ncstore)
!mydist u(*,*,block,*),vt(*,block,*,*,*),res(*,block,*)
!mydist c(*,block,*)

!
! locals
integer i,j,nf,l,lk,im,jm,kp,ip,jp,km, nf5, lkp, lkm
parameter(nf5=5)
! 5 is the number of buffer zones for a domain of dependence of 2.
real  vk(5,maxfluidcmp,ng1,ng2),&
      vij(ng1,ng2,maxfluidcmp), u1(5,maxfluidcmp,ng1,ng2),&
      u2(ng1,ng2,maxfluidcmp), flz(2,ng1,ng2,maxfluidcmp),&
      flx(ng1,ng2,maxfluidcmp),fp,fm
!mydist vk(*,*,*,block),vij(*,block,*),u1(*,*,*,block)
!mydist u2(*,block,*),flz(*,*,block,*),flx(*,block,*)
! to conserve memory, one could equivalence some of the arrays.  But
! I am worried that equivalences break many optimizations and
! parallelizations.
integer lastk
save lastk
data lastk /0/
#ifdef COLD    
      include 'cold.f90'
#endif      

!dec$ alias addr_to_rad, "_OtsAddrToRad"
integer  addr_to_rad,irad,irad1,irad2
external addr_to_rad
integer omp_get_thread_num,cpuid,cpu_get_rad
external omp_get_thread_num,cpuid,cpu_get_rad

!cdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(vk,vij,u1,u2,flz,flx)
           
if (nfluid .ne. nf5) then
   write(*,*) 'dpv: only nfluid=5 currently support, nf=',nfluid
   stop
endif
! first sweep the k dimension:
! zpj: I make c$doacross comments because it interfere with c$omp
!c$doacross local(j,nf,l,i,lk)
!*KAP*parallel region local(j,nf,l,i,lk)
!*KAP*parallel do

!$omp parallel do default(private) firstprivate(k) & 
!$omp shared(vk,vt,c,indxr,indxc,u,u1,indxu)
do j=1,ng2

   do nf=1,nf5
! l covers the domain of dependence
      do l=1,5
         lk = mod(k+l-4+ng3,ng3)+1
         do i=1,ng1
            vk(l,nf,i,j)=vt(i,j,indxr(lk),3,nf)/c(i,j,indxc(lk))
         enddo
      enddo
   enddo
   do l=1,5
   lk = mod(k+l-4+ng3,ng3)+1
      do i=1,ng1
            u1(l,1,i,j)=u(1,i,j,indxu(lk))
            u1(l,2,i,j)=u(2,i,j,indxu(lk))
            u1(l,3,i,j)=u(3,i,j,indxu(lk))
            u1(l,4,i,j)=u(4,i,j,indxu(lk))
            u1(l,5,i,j)=u(5,i,j,indxu(lk))
      enddo
   enddo
#if 0
      if (k .eq. 1) then
         irad1=addr_to_rad(u1(1,1,1,j))
         irad2=addr_to_rad(u(1,1,j,1))
         irad=cpu_get_rad()
         if ((irad1 .ne. irad) .or. (irad2 .ne. irad)) then
                 write(*,*) j,irad1,irad2,irad
         endif
      endif
#endif
enddo
!$omp end parallel do

!*KAP*end parallel region
call wphalfz(flz,u1,vk,nfluid)

kp=mod(k,ng3)+1
km=mod(k+ng3-2,ng3)+1
lk = mod(k+3-4+ng3,ng3)+1
lkm=mod(k+2-4+ng3,ng3)+1
lkp=mod(k+4-4+ng3,ng3)+1
      
!c*$* assert do (concurrent)

! zpj: I make c$doacross comments
!c$doacross local(j,nf,i) shared(ng2,ng1,res,indxc,kp,k,c,km,flz
!c$&  ,vij,vt,indxr,u,u2,indxu)

!*KAP*parallel region local(j,nf,i)
!*KAP*parallel do
! zpj: end
#ifdef COLD
!$omp parallel do default(private) & 
!$omp shared(flz,vt,indxr,indxc,c,res,vij,u2,u,indxu,cold,k,kp,km,lk,lkm,lkp)
#else
!$omp parallel do default(private) &
!$omp shared(flz,vt,indxr,indxc,c,res,vij,u2,u,indxu,k,kp,km,lk,lkm,lkp)
#endif 
do j=1,ng2
#ifdef COLD         
   do i=1,ng1
      if (cold(i,j,k) .and. cold(i,j,kp) ) then
         do nf=1,5
            flz(2,i,j,nf)=(vt(i,j,indxr(lk),3,nf)+&
      vt(i,j,indxr(lkp),3,nf))/(c(i,j,indxc(kp))+c(i,j,indxc(k)))
         enddo
      endif
      if (cold(i,j,km) .and. cold(i,j,k) ) then
         do nf=1,5
            flz(1,i,j,nf)=(vt(i,j,indxr(lk),3,nf)+&
      vt(i,j,indxr(lkm),3,nf))/(c(i,j,indxc(km))+c(i,j,indxc(k)))
         enddo
      endif
   enddo
#endif            
   do nf=1,5
      do i=1,ng1
         res(i,j,nf)=(c(i,j,indxc(kp))+c(i,j,indxc(k)))/2*&
              flz(2,i,j,nf)-(c(i,j,indxc(k))+c(i,j,indxc(km)))/2&
              *flz(1,i,j,nf)
         vij(i,j,nf)=vt(i,j,indxr(k),1,nf)/c(i,j,indxc(k))
         u2(i,j,nf)=u(nf,i,j,indxu(k))
      enddo
   enddo
enddo
!$omp end parallel do
!*KAP*end parallel region
call wphalfx(flx,u2,vij,nfluid)
!c*$* assert do (concurrent)s
! zpj: I make c$doacross comment
!c$doacross local(j,nf,i,ip,im,fp,fm)
#ifdef COLD
!$omp parallel do default(private) & 
!$omp shared(flx,vt,indxr,c,indxc,res,vij,cold,k,lk)
#else
!$omp parallel do default(private) & 
!$omp shared(flx,vt,indxr,c,indxc,res,vij,k,lk)
#endif
do j=1,ng2
#ifdef COLD         
   do i=1,ng1
      ip=mod(i,ng1)+1
      im=mod(i-2+ng1,ng1)+1
      if (cold(i,j,k) .and. cold(ip,j,k) ) then
         do nf=1,5
            flx(i,j,nf)=(vt(i,j,indxr(lk),1,nf)+&
      vt(ip,j,indxr(lk),1,nf))/(c(ip,j,indxc(k))+c(i,j,indxc(k)))
         enddo
      endif
   enddo
#endif            
   do nf=1,nf5
      do i=1,ng1
         ip=mod(i,ng1)+1
         im=mod(i-2+ng1,ng1)+1
         fp=(c(ip,j,indxc(k))+c(i,j,indxc(k)))/2*flx(i,j,nf)
         fm=(c(i,j,indxc(k))+c(im,j,indxc(k)))/2*flx(im,j,nf)
         res(i,j,nf)=res(i,j,nf)+fp-fm
      enddo
   enddo

! lastly the y dimension flux:
   do nf=1,nf5
      do i=1,ng1
         vij(i,j,nf)=vt(i,j,indxr(k),2,nf)/c(i,j,indxc(k))
      enddo
   enddo
enddo
!$omp end parallel do
call wphalfy(flx,u2,vij,nfluid)
#ifdef COLD
!$omp parallel do default(private) & 
!$omp shared(flx,vt,indxr,c,indxc,res,cold,k,lk)
#else
!$omp parallel do default(private) &
!$omp shared(flx,vt,indxr,c,indxc,res,k,lk)
#endif
#ifdef COLD         
!*$* assert do (concurrent)
do j=1,ng2
   jp=mod(j,ng2)+1
   jm=mod(j-2+ng2,ng2)+1
   do i=1,ng1
      if (cold(i,j,k) .and. cold(i,jp,k) ) then
         do nf=1,5
            flx(i,j,nf)=(vt(i,j,indxr(lk),2,nf)+&
      vt(i,jp,indxr(lk),2,nf))/(c(i,jp,indxc(k))+c(i,j,indxc(k)))
         enddo
      endif
   enddo
enddo
!$omp end parallel do
!$omp parallel do default(private) & 
!$omp shared(flx,vt,indxr,c,indxc,res,cold,k,lk)
#endif            
!*$* assert do (concurrent)
do j=1,ng2
   do nf=1,nf5
      do i=1,ng1
         jp=mod(j,ng2)+1
         jm=mod(j-2+ng2,ng2)+1
         res(i,j,nf)=res(i,j,nf)&
              +(c(i,jp,indxc(k))+c(i,j,indxc(k)))/2*flx(i,j,nf)&
              -(c(i,j,indxc(k))+c(i,jm,indxc(k)))/2*flx(i,jm,nf)
      enddo
   enddo
enddo
!$omp end parallel do

if (k .gt. lastk) then
!         write(*,'(5f15.9)') (flz(1,1,2,j),j=1,5)
   lastk=k
else
!          write(*,'(5f15.9)') (vk(1,1,5,j),j=1,5)
endif
return
end

!**********************************************************************
subroutine wphalfx(flx,u,v,nfluid)
! the result flx is one half cell offset to the right.  So
! flx(1)<->u(1.5)
implicit none
#include 'relaxgl.f90'
integer nfluid
real u(ng1,ng2,nfluid),v(ng1,ng2,nfluid),flx(ng1,ng2,nfluid)
!mydist u(*,block,*),v(*,block,*),flx(*,block,*)
! locals
integer i,j,ip,im,ipp,nf
real a,b,ap,bp,da,dam,db,dbm, minmod,mm1,mm2


! minmod is a real misnomer, since it is only minmod if wslopelim=1 or wbee=0.
!
mm1(a,b)=min(wslopelim*a,wbee*b+(1-wbee)*a)
mm2(a,b)=mm1(min(a,b),max(a,b))
minmod(a,b)=(sign(1.,a)+sign(1.,b))*mm2(abs(a),abs(b))/2


!$omp parallel default(private) shared(u,v,flx,nfluid)
!fpp$ skip
!c*$* assert do(serial)      
do nf=1,nfluid
!$dir prefer_parallel_ext            
!c$doacross local(i,j)
!*$* assert do (concurrent)
!$omp do
   do j=1,ng2
      do i=1,ng1
         ip=mod(i,ng1)+1
         ipp=mod(ip,ng1)+1
         im=mod(i+ng1-2,ng1)+1
         a=(u(i,j,nf)+v(i,j,nf))/2
         b=(u(ip,j,nf)-v(ip,j,nf))/2
         da=(u(ip,j,nf)+v(ip,j,nf))/2-a
         dam=a-(u(im,j,nf)+v(im,j,nf))/2
         ap=a+minmod(dam,da)/2
         db=(u(ipp,j,nf)-v(ipp,j,nf))/2-b
         dbm=b-(u(i,j,nf)-v(i,j,nf))/2
         bp=b-minmod(db,dbm)/2
         flx(i,j,nf)=ap-bp
      enddo
   enddo
!$omp enddo nowait
enddo
!$omp end parallel
return
end

!**********************************************************************
subroutine wphalfy(flx,u,v,nfluid)
! the result flx is one half cell offset to the right.  So
! flx(1)<->u(1.5)
implicit none
#include 'relaxgl.f90'
integer nfluid
real u(ng1,ng2,nfluid),v(ng1,ng2,nfluid),flx(ng1,ng2,nfluid)
!mydist u(*,block,*),v(*,block,*),flx(*,block,*)

! locals
integer i,j,jp,jm,jpp,nf
real a,b,ap,bp,da,dam,db,dbm,minmod,mm1,mm2

mm1(a,b)=min(wslopelim*a,wbee*b+(1-wbee)*a)
mm2(a,b)=mm1(min(a,b),max(a,b))
minmod(a,b)=(sign(1.,a)+sign(1.,b))*mm2(abs(a),abs(b))/2


!$omp parallel default(private) shared(nfluid,flx,u,v)
!fpp$ skip
!c*$* assert do(serial)      
do nf=1,nfluid
!$dir prefer_parallel_ext            
!c$doacross local(i,j)
!*$* assert do (concurrent)
!$omp do
   do j=1,ng2
      jp=mod(j,ng2)+1
      jpp=mod(jp,ng2)+1
      jm=mod(j+ng2-2,ng2)+1
      do i=1,ng1
         a=(u(i,j,nf)+v(i,j,nf))/2
         b=(u(i,jp,nf)-v(i,jp,nf))/2
         da=(u(i,jp,nf)+v(i,jp,nf))/2-a
         dam=a-(u(i,jm,nf)+v(i,jm,nf))/2
         ap=a+minmod(dam,da)/2
         db=(u(i,jpp,nf)-v(i,jpp,nf))/2-b
         dbm=b-(u(i,j,nf)-v(i,j,nf))/2
         bp=b-minmod(db,dbm)/2
         flx(i,j,nf)=ap-bp
      enddo
   enddo
!$omp enddo nowait
enddo
!$omp end parallel
return
end



!**********************************************************************
subroutine wphalfz(flz,u,v,nfluid)
! the result flx is one half cell offset to the right.  So
! flx(1)<->u(1.5)
implicit none
#include 'relaxgl.f90'
integer nfluid
real u(5,nfluid,ng1,ng2),v(5,nfluid,ng1,ng2),flz(2,ng1,ng2,nfluid)
!mydist u(*,*,*,block),v(*,*,*,block),flz(*,*,block,*)

! locals
integer i,j,nf,k1,k,km,kp,kpp
real a,b,ap,bp,da,dam,db,dbm,minmod,mm1,mm2

mm1(a,b)=min(wslopelim*a,wbee*b+(1-wbee)*a)
mm2(a,b)=mm1(min(a,b),max(a,b))
minmod(a,b)=(sign(1.,a)+sign(1.,b))*mm2(abs(a),abs(b))/2

!$omp parallel do default(private) shared(u,v,flz)
!*$* assert do (concurrent)
do j=1,ng2
   do i=1,ng1
      nf=1
      k1=1
      k=k1+1
      km=k1
      kp=k1+2
      kpp=k1+3
            a=(u(k,nf,i,j)+v(k,nf,i,j))/2
            b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
            da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
            dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
            ap=a+minmod(dam,da)/2
            db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
            dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
            bp=b-minmod(db,dbm)/2
            flz(k1,i,j,nf)=ap-bp
      k1=2
      k=k1+1
      km=k1
      kp=k1+2
      kpp=k1+3
            a=(u(k,nf,i,j)+v(k,nf,i,j))/2
            b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
            da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
            dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
            ap=a+minmod(dam,da)/2
            db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
            dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
            bp=b-minmod(db,dbm)/2
            flz(k1,i,j,nf)=ap-bp
      nf=2
      k1=1
      k=k1+1
      km=k1
      kp=k1+2
      kpp=k1+3
            a=(u(k,nf,i,j)+v(k,nf,i,j))/2
            b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
            da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
            dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
            ap=a+minmod(dam,da)/2
            db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
            dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
            bp=b-minmod(db,dbm)/2
            flz(k1,i,j,nf)=ap-bp
      k1=2
      k=k1+1
      km=k1
      kp=k1+2
      kpp=k1+3
            a=(u(k,nf,i,j)+v(k,nf,i,j))/2
            b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
            da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
            dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
            ap=a+minmod(dam,da)/2
            db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
            dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
            bp=b-minmod(db,dbm)/2
            flz(k1,i,j,nf)=ap-bp
      nf=3
      k1=1
      k=k1+1
      km=k1
      kp=k1+2
      kpp=k1+3
            a=(u(k,nf,i,j)+v(k,nf,i,j))/2
            b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
            da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
            dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
            ap=a+minmod(dam,da)/2
            db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
            dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
            bp=b-minmod(db,dbm)/2
            flz(k1,i,j,nf)=ap-bp
      k1=2
      k=k1+1
      km=k1
      kp=k1+2
      kpp=k1+3
            a=(u(k,nf,i,j)+v(k,nf,i,j))/2
            b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
            da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
            dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
            ap=a+minmod(dam,da)/2
            db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
            dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
            bp=b-minmod(db,dbm)/2
            flz(k1,i,j,nf)=ap-bp
      nf=4
      k1=1
      k=k1+1
      km=k1
      kp=k1+2
      kpp=k1+3
            a=(u(k,nf,i,j)+v(k,nf,i,j))/2
            b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
            da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
            dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
            ap=a+minmod(dam,da)/2
            db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
            dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
            bp=b-minmod(db,dbm)/2
            flz(k1,i,j,nf)=ap-bp
      k1=2
      k=k1+1
      km=k1
      kp=k1+2
      kpp=k1+3
            a=(u(k,nf,i,j)+v(k,nf,i,j))/2
            b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
            da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
            dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
            ap=a+minmod(dam,da)/2
            db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
            dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
            bp=b-minmod(db,dbm)/2
            flz(k1,i,j,nf)=ap-bp
      nf=5
      k1=1
      k=k1+1
      km=k1
      kp=k1+2
      kpp=k1+3
            a=(u(k,nf,i,j)+v(k,nf,i,j))/2
            b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
            da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
            dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
            ap=a+minmod(dam,da)/2
            db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
            dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
            bp=b-minmod(db,dbm)/2
            flz(k1,i,j,nf)=ap-bp
      k1=2
      k=k1+1
      km=k1
      kp=k1+2
      kpp=k1+3
            a=(u(k,nf,i,j)+v(k,nf,i,j))/2
            b=(u(kp,nf,i,j)-v(kp,nf,i,j))/2
            da=(u(kp,nf,i,j)+v(kp,nf,i,j))/2-a
            dam=a-(u(km,nf,i,j)+v(km,nf,i,j))/2
            ap=a+minmod(dam,da)/2
            db=(u(kpp,nf,i,j)-v(kpp,nf,i,j))/2-b
            dbm=b-(u(k,nf,i,j)-v(k,nf,i,j))/2
            bp=b-minmod(db,dbm)/2
            flz(k1,i,j,nf)=ap-bp
   enddo
enddo
return
end


!**********************************************************************
subroutine stepu(unew, uold, res, dt, indxunew,  ngul3, nfluidcmp, k)
implicit none
#include 'relaxgl.f90'
integer k, indxunew(ng3), ngul3, nfluidcmp
real dt, unew(nfluidcmp,ng1,ng2,ngul3),& 
     uold(nfluidcmp,ng1,ng2,ng3), res(ng1,ng2,nfluidcmp)
!mydist unew(*,*,block,*),uold(*,*,block,*)
!mydist res(*,block,*)
!
! locals
integer i,j,nf

!fpp$ skip
!c*$* assert do(serial)  
!cdir loopcnt=5
!$omp parallel private(i,j) firstprivate(nfluidcmp) &
!$omp shared(unew,indxunew,uold,k,res,dt)   
do nf=1,nfluidcmp
!$dir prefer_parallel_ext            
!c$doacross local(i,j)
!*$* assert do (concurrent)
!$omp do
   do j=1,ng2
!cdir vector
      do i=1,ng1
         unew(nf,i,j,indxunew(k))=uold(nf,i,j,k)-res(i,j,nf)*dt
      enddo
   enddo
!$omp enddo nowait
enddo
!$omp end parallel
return
endsubroutine stepu

!**********************************************************************
subroutine stepuu( u, res, dt, nfluidcmp, k)
!
implicit none
#include 'relaxgl.f90'
integer k, nfluidcmp
real dt, u(nfluidcmp,ng1,ng2,ng3), res(ng1,ng2,nfluidcmp)
!mydist u(*,*,block,*), res(*,block,*)
! locals
integer i,j,nf

!fpp$ skip
!c*$* assert do(serial)
!$omp parallel default(private) shared(nfluidcmp,u,res,k,dt)
!cdir loopcnt=5
do nf=1,nfluidcmp
!$dir prefer_parallel_ext
!*$* assert do (concurrent)
!$omp do
    do j=1,ng2
!cdir vector
      do i=1,ng1
         u(nf,i,j,k)=u(nf,i,j,k)-res(i,j,nf)*dt
      enddo
   enddo
!$omp enddo nowait
enddo
!$omp end parallel
return
endsubroutine stepuu
      
!**********************************************************************
subroutine calcbaj(baj,def,indx,k)
implicit none
#include 'relaxgl.f90'
#ifdef GMETRIC
include 'gmetric.fi'
#endif
integer k,indx(ng3)
real baj(7,ng1,ng2,nstore),def(ng1,ng2,ng3)
!mydist baj(*,*,block,*),def(*,block,*)
!
! locals
integer i,j,ip,jp,im,jm,kp,km
real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22,&
     a23, a33, det, b11,b12,b13,b22,b23,b33


!$omp parallel do default(private) shared(baj,k,indx,def)
!$dir prefer_parallel_ext            
!c$doacross local(i,j)
do j=1,ng2
   do i=1,ng1
#ifdef GMETRIC
      baj(1,i,j,indx(k))=gbaj(1,i,j,k)
      baj(2,i,j,indx(k))=gbaj(2,i,j,k)
      baj(3,i,j,indx(k))=gbaj(3,i,j,k)
      baj(4,i,j,indx(k))=gbaj(4,i,j,k)
      baj(5,i,j,indx(k))=gbaj(5,i,j,k)
      baj(6,i,j,indx(k))=gbaj(6,i,j,k)
      baj(7,i,j,indx(k))=gbaj(7,i,j,k)
#else
      ip=mod(i,ng1)+1
      im=mod(i+ng1-2,ng1)+1
      jp=mod(j,ng2)+1
      jm=mod(j+ng2-2,ng2)+1
      kp=mod(k,ng3)+1
      km=mod(k+ng3-2,ng3)+1
      phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))
      phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))
      phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))
      phixy=(def(ip,jp,k)-def(im,jp,k)-def(ip,jm,k)+def(im,jm,k))/4
      phiyz=(def(i,jp,kp)-def(i,jp,km)-def(i,jm,kp)+def(i,jm,km))/4
      phixz=(def(ip,j,kp)-def(im,j,kp)-def(ip,j,km)+def(im,j,km))/4
      a11=1+phixx
      a12=phixy
      a13=phixz
      a22=1+phiyy
      a23=phiyz
      a33=1+phizz
      det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2-a12**2*a33
      baj(7,i,j,indx(k))=det
      b11=(a22*a33-a23**2)/det
      b12=(a13*a23-a12*a33)/det
      b13=(a12*a23-a13*a22)/det
      b22=(a11*a33-a13**2)/det
      b23=(a12*a13-a11*a23)/det
      b33=(a11*a22-a12**2)/det
      baj(1,i,j,indx(k))=b11
      baj(2,i,j,indx(k))=b12
      baj(3,i,j,indx(k))=b13
      baj(4,i,j,indx(k))=b22
      baj(5,i,j,indx(k))=b23
      baj(6,i,j,indx(k))=b33
#endif
   enddo
enddo
return
endsubroutine calcbaj


!**********************************************************************
subroutine flux(vt,defp,u,baj,indxb,indxu,nfluidcmp,ngu3,k)
implicit none
#include 'relaxgl.f90'
integer k, nfluidcmp, indxb(ng3),indxu(ng3),ngu3
real vt(ng1,ng2,nstore,3,nfluidcmp),defp(ng1,ng2,ng3),&
     baj(7,ng1,ng2,nstore),u(nfluidcmp,ng1,ng2,ngu3),&
     dascale
external dascale
!mydist vt(*,block,*,*,*),defp(*,block,*)
!mydist baj(*,*,block,*),u(*,*,block,*)
!
!     locals
real rootg(ng1,ng2),vt1(3,maxfluidcmp,ng1,ng2),&
     ul(maxfluidcmp,ng1,ng2),rg,pdxg,pdyg,pdzg,&
     b11,b12,b13,b22,b23,b33,vtt1,vtt2,vtt3,&
     pdcache(29,ng1,ng2),t
!mydist rootg(*,block),vt1(*,*,*,block),ul(*,*,block)
!mydist pdcache(*,*,block)
integer i,j,ndim,nf,ip,jp,kp,im,jm,km,l, ik
integer lastk
save lastk
data lastk /0/

ik=indxb(k)
! copy the array into a temporary for efficiency and portability
! across distributed memory machines and f77 compilers.
! 7.2 f77
! zpj: I made c$doacros comments
!c$doacross local(i,nf)

!$omp parallel do default(private) & 
!$omp shared(rootg,baj,ul,u,indxu,ik,k,nfluidcmp)
do j=1,ng2
!cdir vector
   do i=1,ng1
      rootg(i,j)=baj(7,i,j,ik)
!cdir loopcnt=5 nounroll assert(nfluidcmp=5)
      do nf=1,nfluidcmp
         ul(nf,i,j)=u(nf,i,j,indxu(k))
      enddo
   enddo
enddo 
call fluxeuler(vt1,ul,rootg,nfluidcmp,k)
      
! 7.2 f77
! zpj: make  c$doacros comments
!c$doacross local(i,ip,im,jp,jm,kp,km,rg,b11,b12,b13,b22,b23,b33,pdxg
!c$&     ,pdyg,pdzg,vtt1,vtt2,vtt3,l)
!$omp parallel do default(private) & 
!$omp shared(rootg,defp,baj,vt1,u,indxu,vt,nfluidcmp,ik,k)
!cdir nodep
do j=1,ng2
!cdir vector, nodep
   do i=1,ng1
      ip=mod(i,ng1)+1
      im=mod(i+ng1-2,ng1)+1
      jp=mod(j,ng2)+1
      jm=mod(j+ng2-2,ng2)+1
      kp=mod(k,ng3)+1
      km=mod(k+ng3-2,ng3)+1
      rg=rootg(i,j)
      pdxg=(defp(ip,j,k)-defp(im,j,k))/2/rg
      pdyg=(defp(i,jp,k)-defp(i,jm,k))/2/rg
      pdzg=(defp(i,j,kp)-defp(i,j,km))/2/rg
      b11=baj(1,i,j,ik)*rg
      b12=baj(2,i,j,ik)*rg
      b13=baj(3,i,j,ik)*rg
      b22=baj(4,i,j,ik)*rg
      b23=baj(5,i,j,ik)*rg
      b33=baj(6,i,j,ik)*rg
!cdir loopcnt=5
      do l=1,nfluidcmp
         vtt1=vt1(1,l,i,j)-u(l,i,j,indxu(k))*pdxg
         vtt2=vt1(2,l,i,j)-u(l,i,j,indxu(k))*pdyg
         vtt3=vt1(3,l,i,j)-u(l,i,j,indxu(k))*pdzg
         vt(i,j,ik,1,l)=b11*vtt1+b12*vtt2+b13*vtt3
         vt(i,j,ik,2,l)=b12*vtt1+b22*vtt2+b23*vtt3
         vt(i,j,ik,3,l)=b13*vtt1+b23*vtt2+b33*vtt3
      enddo
   enddo
enddo
if (k .gt. lastk) then
!         write(*,'(3f20.9)') (vt1(1,1,j,3),j=1,3)
   lastk=k
endif
return
endsubroutine flux

!**********************************************************************
subroutine fluxeuler(vt,u,rootg,nfluidcomp,k)
implicit none
#include 'relaxgl.f90'
integer nfluidcomp, k
real vt(3,nfluidcomp,ng1,ng2), u(nfluidcomp,ng1,ng2),rootg(ng1,ng2)
!mydist vt(*,*,*,block),u(*,*,block),rootg(*,block)
!
! locals
integer i,j
real rho,mx,my,mz,vx,vy,vz,engy,kinetic,pressure,econv, gamma, rg
parameter(gamma=5./3.)
#ifdef COLD
      include 'cold.f90'
#endif      


if (nfluidcomp .ne. 5) then
   write(*,*) 'fluxeuler: invalid argument: nfluid=',nfluidcomp
   stop
endif

!$dir prefer_parallel_ext            
!c$doacross local(i,j)
#ifdef COLD
!$omp parallel do default(private) shared(rootg,u,vt,k,cold)
#else
!$omp parallel do default(private) shared(rootg,u,vt,k)
#endif
do j=1,ng2
   do i=1,ng1
      rg=rootg(i,j)
      rho=u(1,i,j)/rg
      mx=u(2,i,j)/rg
      my=u(3,i,j)/rg
      mz=u(4,i,j)/rg
      vx=mx/rho
      vy=my/rho
      vz=mz/rho
      engy=u(5,i,j)/rg
      kinetic=rho*(vx**2+vy**2+vz**2)/2
      pressure=(gamma-1)*(engy-kinetic)
      pressure=max(0.,pressure)
#ifdef  COLD
      if (cold(i,j,k)) pressure=0
#endif            
      econv=engy+pressure
      vt(1,1,i,j)=mx
      vt(2,1,i,j)=my
      vt(3,1,i,j)=mz
      vt(1,2,i,j)=mx*vx+pressure
      vt(2,2,i,j)=mx*vy
      vt(3,2,i,j)=mx*vz
      vt(1,3,i,j)=my*vx
      vt(2,3,i,j)=my*vy+pressure
      vt(3,3,i,j)=my*vz
      vt(1,4,i,j)=mz*vx
      vt(2,4,i,j)=mz*vy
      vt(3,4,i,j)=mz*vz+pressure
      vt(1,5,i,j)=econv*vx
      vt(2,5,i,j)=econv*vy
      vt(3,5,i,j)=econv*vz
   enddo
     enddo
return
end         


!**********************************************************************
subroutine cflkeuler(c,u,baj,defp,indxc,indx,k)
!
! to save on memory, we calculate the metric from scratch.
!
implicit none
#include 'relaxgl.f90'
integer k,indx(ng3),indxc(ng3)
real c(ng1,ng2,ncstore),u(5,ng1,ng2,ng3),baj(7,ng1,ng2,nstore),defp(ng1,ng2,ng3)
!mydist c(*,block,*),u(*,*,block,*),baj(*,*,block,*)
!mydist defp(*,block,*)
! locals
integer i,j,ip,im,jp,jm,kp,km
real det,pdx,pdy,pdz,vx,vy,vz,rho,b11,b12,b13,b22,b23,b33,&
     kin,ed,pressure,cs,bnorm,tcfl, blim1, blim2, blim3,&
     vv1,vv2,vv3, v1, v2, v3
real gamma
parameter(gamma=5./3.)

!$omp parallel do default(private) & 
!$omp shared(defp,baj,u,c,indxc,indx,k)
!$dir prefer_parallel_ext         
!c$doacross local(i,j)
do j=1,ng2
   do i=1,ng1
      ip=mod(i,ng1)+1
      im=mod(i+ng1-2,ng1)+1
      jp=mod(j,ng2)+1
      jm=mod(j+ng2-2,ng2)+1
      kp=mod(k,ng3)+1
      km=mod(k+ng3-2,ng3)+1
      pdx=(defp(ip,j,k)-defp(im,j,k))/2
      pdy=(defp(i,jp,k)-defp(i,jm,k))/2
      pdz=(defp(i,j,kp)-defp(i,j,km))/2
      b11=baj(1,i,j,indx(k))
      b12=baj(2,i,j,indx(k))
      b13=baj(3,i,j,indx(k))
      b22=baj(4,i,j,indx(k))
      b23=baj(5,i,j,indx(k))
      b33=baj(6,i,j,indx(k))
      det=baj(7,i,j,indx(k))
      vx=u(2,i,j,k)/u(1,i,j,k)
      vy=u(3,i,j,k)/u(1,i,j,k)
      vz=u(4,i,j,k)/u(1,i,j,k)
      rho=u(1,i,j,k)/det
      kin=rho*(vx**2+vy**2+vz**2)/2
      ed=u(5,i,j,k)/det-kin
      pressure=(gamma-1)*ed
      cs=sqrt(gamma*abs(pressure/rho))
! one should really calculate the largest eigenvalue of the
! dreibein, but I''m too lazy, so will use the infinity norm
! as a bound.
!            bnorm=max(abs(b11)+abs(b12)+abs(b13),abs(b12)
!     &           +abs(b22)+abs(b23),abs(b13),abs(b23),abs(b33))
!            v2=sqrt((vx-pdx)**2+(vy-pdy)**2+(vz-pdz)**2)
!            tcfl=bnorm*(v2+cs)

! this interesting formula is apparently an exact result, as
! listed in Yee and confirmed by gnedin.
      v1=vx-pdx
      v2=vy-pdy
      v3=vz-pdz
      vv1=b11*v1+b12*v2+b13*v3
      vv2=b12*v1+b22*v2+b23*v3
      vv3=b13*v1+b23*v2+b33*v3
      blim1=cs*sqrt(b11**2+b12**2+b13**2)+abs(vv1)
      blim2=cs*sqrt(b12**2+b22**2+b23**2)+abs(vv2)
      blim3=cs*sqrt(b13**2+b23**2+b33**2)+abs(vv3)
      tcfl=max(blim1,blim2,blim3)
      c(i,j,indxc(k))=tcfl
   enddo
enddo
return
endsubroutine 

