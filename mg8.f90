!-*- Fortran -*-	

subroutine mg8(u,rhs,defpot,ubig,dx,nu,iopt)
implicit none
integer iopt,nu
real u(8,8,8),rhs(8,8,8),defpot(8,8,8),ubig(nu,8,8,8),dx
!dir$ shared *u(:block,:block,:),*defpot(:block,:block,:)
!dir$ shared *rhs(:block,:block,:)      
! locals
integer isolver,imgalg,irelaxopt,nrelax,i
#define N2 8/2
real uhalf(N2,N2,N2),rhshalf(N2,N2,N2),defhalf(N2,N2,N2),usmall(1,N2,N2,N2)
!dir$ shared uhalf(:block,:block,:),defhalf(:block,:block,:)
!dir$ shared rhshalf(:block,:block,:),usmall(:,:block,:block,:)
real relaxerr,r1

isolver=mod(iopt,8)

#ifdef _ALPHA
! first touch placement of the subarrays
!cdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(uhalf,defhalf,rhshalf,usmall)
call szero(uhalf,N2,N2,N2)
call szero(defhalf,N2,N2,N2)
call szero(rhshalf,N2,N2,N2)
call szero(usmall,N2,N2,N2)
#endif
nrelax=4*NG/8
irelaxopt=isolver*8+3
! calculate residual on half grid      
call relax(u,rhshalf,rhs,defpot,ubig,r1,dx,8,8,8,nu,nrelax,irelaxopt) 

! goto 59
call restrict(defhalf,defpot,8,8,8)
call restrictu(usmall,ubig,8,8,8,nu)
call szero(uhalf,N2,N2,N2)

call mg4(uhalf,rhshalf,defhalf,usmall,dx*2,1,iopt)

call inject(u,uhalf,8,8,8)
59    continue
irelaxopt=isolver*8+2
call relax(u,rhshalf,rhs,defpot,ubig,relaxerr,dx,8,8,8,nu,nrelax,irelaxopt) 
if (8 .eq. NG) then
   write(*,*) 'relaxerr=',relaxerr
endif
if (8 .ne. NG .and. relaxerr .gt. r1 ) then
  write(*,73) 8,r1,relaxerr
73      format(' mg: relaxation diverged, level=',i3,' r1,2=',2e12.4)
! prolong resets u: forget about the current level relaxation!
!  write(*,*) 'call prolong 8'
!  call prolong(u,uhalf,8,8,8)
  call szero(u,8,8,8)
else
  call zerosum(u,8,8,8)
endif
return
endsubroutine mg8
