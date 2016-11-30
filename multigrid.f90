! -*- fortran -*-
#include "dimen.f90"

#ifdef TESTMG
call testmg
end
#endif
      
!      
! f77 static memory version:
! multigrid solver, based on gauss-seidel relaxation.
!
! written November 11, 1995 by Ue-Li Pen, upen@cfa.harvard.edu
!
! as before, iopt=2 means potential flow solver, iopt=1 is gravity
subroutine multigrid(u,rhs,defpot,ubig,dx,nx,ny,nz,nu,nred,iopt)
implicit none
integer nx,ny,nz,iopt,nred,nu
real u(nx,ny,nz),defpot(nx,ny,nz),rhs(nx,ny,nz),dx,ubig(nu,nx,ny,nz),asum
!dir$ shared *u(:block,:block,:),*defpot(:block,:block,:)
!dir$ shared *rhs(:block,:block,:), *ubig(:,:block,:block,:)

integer i5
! we will call specific routines for each power of 2 dimension.



if (nx .ne. ny .or. ny .ne. nz) then
   write(*,*) 'multigrid: only perfect cubes supported'
   stop
endif
if (nred .ne. 1) then
   write(*,*) 'multigrid: only nred=1 supported, nred=',nred
   stop
endif

call rgzerosum(rhs,defpot,dx,nx,ny,nz)
!      if (mod(iopt,4) .eq. 1) call rgmult(u,defpot,nx,ny,nz)
! we have to set the mean to zero here.  In actuality, if u was
! unchanged since the last call to multigrid, this is superfluous.
call usum(asum,u,nx,ny,nz)
if (asum .gt. 1.e-2) then
   write(*,*) 'multigrid: u is not zero summed, calling zerosum'
   write(*,*) 'sum(u)=',asum
   call zerosum(u,nx,ny,nz)
endif

do i5=1,1
if (nx .eq. 4) then
   call mg4(u,rhs,defpot,ubig,dx,nu,iopt)
#if NG>=8
else if (nx .eq. 8) then
   call mg8(u,rhs,defpot,ubig,dx,nu,iopt)
#endif
#if NG>=16         
else if (nx .eq. 16) then
   call mg16(u,rhs,defpot,ubig,dx,nu,iopt)
#endif
#if NG>=32         
else if (nx .eq. 32) then
   call mg32(u,rhs,defpot,ubig,dx,nu,iopt)
#endif
#if NG>=64
else if (nx .eq. 64) then
   call mg64(u,rhs,defpot,ubig,dx,nu,iopt)
#endif
#if NG>=128         
else if (nx .eq. 128) then
   call mg128(u,rhs,defpot,ubig,dx,nu,iopt)
#endif
#if NG>=256        
else if (nx .eq. 256) then
   call mg256(u,rhs,defpot,ubig,dx,nu,iopt)
#endif
#if NG>=512        
else if (nx .eq. 512) then
   call mg512(u,rhs,defpot,ubig,dx,nu,iopt)
#endif         
#if NG>=1024        
else if (nx .eq. 1024) then
   call mg1024(u,rhs,defpot,ubig,dx,nu,iopt)
#endif         
else
   write(*,*) 'multigrid: nx=',nx,'  not supported'
   stop
endif
enddo
!      if (mod(iopt,4) .eq. 1) call rgdiv(u,defpot,nx,ny,nz)
return
end


      
subroutine mg4(u,rhs,defpot,ubig,dx,nu,iopt)
implicit none
integer ng,nu
parameter(ng=4)
real u(ng,ng,ng),rhs(ng,ng,ng),defpot(ng,ng,ng),dx,ubig(nu,ng,ng,ng)
!dir$ shared *u(:block,:block,:),*defpot(:block,:block,:)
!dir$ shared *rhs(:block,:block,:), *ubig(:,:block,:block,:)
integer iopt
! locals
integer isolver,imgalg,irelaxopt,nrelax
real rhshalf(ng/2,ng/2,ng/2),relaxerr,r1
!dir$ shared rhshalf(:block,:block,:)
      
isolver=mod(iopt,8)
nrelax=32
irelaxopt=isolver*8+2
call relax(u,rhshalf,rhs,defpot,ubig,r1,dx,ng,ng,ng,nu,nrelax,irelaxopt) 
nrelax=32
call relax(u,rhshalf,rhs,defpot,ubig,relaxerr,dx,ng,ng,ng,nu,nrelax,irelaxopt) 
if (ng .ne. NG .and. relaxerr .gt. r1*1.001 .and. r1 .gt. (1.E-5)) then
  write(*,73) ng,r1,relaxerr
73      format(' mg: relaxation diverged, level=',i3,' r1,2=',2e12.4)
  call szero(u,ng,ng,ng)
else
  call zerosum(u,ng,ng,ng)
endif
return
endsubroutine mg4

#ifdef TESTMG

subroutine testmg
real u(NG,NG,NG),rhs(NG,NG,NG),def(NG,NG,NG),t(NG),t2(NG),ubig(1,NG,NG,NG),t3(NG),t4(NG)
!dir$ shared u(:block,:block,:),def(:block,:block,:)
!dir$ shared rhs(:block,:block,:),ubig(:,:block,:block,:)
! test solver in linear 1-D:
! set ubig to be a sine wave:

def=0
write(*,*) 'testing 3D multigrid, NG=',NG
u=0
pi=4*atan(1.)
ubig(1,:,1,1)=1-(/ (0.5*sin(2*pi*i/NG) ,i=1,NG ) /)
ubig(1,:,:,:)=spread(spread(ubig(1,:,1,1),1,NG),1,NG)
rhs=sqrt(ubig(1,:,:,:))
def=cshift(ubig(1,:,:,:),NG/3)-1
!      rhs(NG/2,:,:)=NG
!      write(*,*) 'sum(rhs)=',sum(rhs)
rhs=rhs-sum(rhs)/NG**3
write(*,*) 'sum(rhs)=',sum(rhs)

do ii=1,10
!dir$ master
!      call secondr(atime)
!dir$ end master
call multigrid(u,rhs,def,ubig,1.,NG,NG,NG,1,1,2)

!dir$ master
!      call secondr(btime)
   write(*,*) 'elapsed wall clock time=',btime-atime
!dir$ end master

t=u(1,1,:)
t4=rhs(1,1,:)
t3=ubig(1,1,1,:)
write(*,*)t
t2=def(1,1,:)
t3=t3/(1-2*t2+cshift(t2,1)+cshift(t2,-1))
s=0
do i=1,NG
   im=mod(i+NG-2,NG)+1
   ip=mod(i,NG)+1
   t2(i)=((t(i)-t(im))*(t3(im)+t3(i))/2-(t(ip)-t(i))*(t3(ip)+t3(i))/2+t4(i))**2
   s=s+t2(i)
enddo
!dir$ master  
write(*,*) sqrt(s/NG)
!dir$ end master      
enddo
return
endsubroutine testmg

#endif
      
#include "mg8.f90"
#include "mg16.f90"
#include "mg32.f90"
#include "mg64.f90"
#include "mg128.f90"
#include "mg256.f90"
#include "mg512.f90"
#include "mg1024.f90"
