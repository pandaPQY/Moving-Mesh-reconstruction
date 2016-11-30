program main
use omp_lib
implicit none
#include "dimen.f90"

character(256),parameter::workdir='/scratch2/p/pen/qiaoyin/datafromsimulation/cafproject/MM_code/'
character(8),parameter::zstring='0.000'
character(64),parameter::inf='den'
character(64),parameter::idstring='_2'
character(8),parameter::prefix='.bin'
character(10),parameter::fnu='u512/'
character(10),parameter::fndef='def512/'
character(10),parameter::fndelta='delta512/'
character(10),parameter::fndisp='disp512/'
integer(4),parameter::niter=300

integer(4),parameter::ng1=NG
integer(4),parameter::ng2=NG
integer(4),parameter::ng3=NG
integer(4),parameter::nfluid=5
real(4)::defp(ng1,ng2,ng3), def(ng1,ng2,ng3),u(nfluid,ng1,ng2,ng3),tmp(ng1,ng2,ng3),tmp2(ng1+2,ng2,ng3)
real(4)::grad(ng1,ng2,ng3,3)
real(4)::delta(ng1,ng2,ng3)
real(4)::umean
integer(4) i,j,k
character(4)::gstring
character(512)::filename
logical(4)::look=.false.
character(40)::fmtstring='(   f12.5)'
! output switcher
logical(4)::flag_def=.true.
logical(4)::flag_gradient=.true.
logical(4)::flag_density=.true.
logical(4)::flag_check=.true.  ! write down checkpoint after last iteration
! restart
logical(4)::restart=.true.
integer(4)::restartstep=150
character(10)::rstring_u='restartu'
character(10)::rstring_d='restartdef'
character(4)::sstring
! misc
integer(4)::iter,begin
real(8)::ompt1,ompt2
real(4)::t1,t2
real(4)::rdummy
include 'globalpa.f90'

! grid issues
write(fmtstring(2:4),'(i3)') ng1
write(gstring,'(i4.4)') NG
write(*,*) 'grid:',NG

write(*,*) ' '
write(*,*) 'BEGIN'

compressmax=10
u=0.
! input density field with mean of 1
call cpu_time(t1)
!filename=trim(workdir)//trim(zstring)//trim(inf)//gstring//trim(idstring)//trim(prefix)
filename=trim(workdir)//trim(zstring)//trim(inf)//gstring//trim(idstring)//trim(prefix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='old',access='stream')
read(31) u(1,:,:,:)
close(31)
call cpu_time(t2)
write(*,*) 'time consumed for reading:',t2-t1,'seconds'

write(*,*) ' '
write(*,*) 'u:'
if (look) write(*,fmtstring) u(1,:,1,1)
umean=sum(real(u(1,:,:,:),8))/float(ng1)/float(ng2)/float(ng3)
!call calcmean(u(1,:,:,:),umean,ng1,ng2,ng3,'u')
if (abs(umean-1.).gt.1.e-5) then
  write(*,*) 'umean should be zero, program ended'
  stop
endif
write(*,*) 'u mean=',real(sum(real(u(1,:,:,:),8))/float(ng1)/float(ng2)/float(ng3))
write(*,*) 'u sigma=',real(sqrt(sum(real((u(1,:,:,:)-1.)**2,8))/float(ng1)/float(ng2)/float(ng3)))

write(*,*) ' '
write(*,*) 'BEGIN MULTIGRID'
def=0.
begin=1
! restart
if (restart) then
   begin=restartstep+1
   call readrestart
endif
do iter=begin,niter
  call cpu_time(t1)
  ompt1=omp_get_wtime()
  write(*,*) ' '
  write(*,*) 'iter=',iter
  call calcdefp(defp,tmp,tmp2,def,u,1.,1.,5)
  write(*,*) 'CALCDEFP'
  write(*,*) 'defp sigma=',real(sqrt(sum(real(defp**2,8))/float(ng1)/float(ng2)/float(ng3)))
  def=def+defp
  write(*,*) 'DEF+DEFP'
  write(*,*) 'def sigma=',real(sqrt(sum(real(def**2,8))/float(ng1)/float(ng2)/float(ng3)))
  call relaxing(u,def,defp,1.,1.,5)
  write(*,*) 'RELAXING'
  write(*,*) 'u sigma=',real(sqrt(sum(real((u(1,:,:,:)-1.)**2,8))/float(ng1)/float(ng2)/float(ng3)))
  write(*,*) 'u min max=',minval(u(1,:,:,:)),maxval(u(1,:,:,:))
  write(*,*) 'u < 0 :',count(u(1,:,:,:).lt.0.)
  call cpu_time(t2)
  ompt2=omp_get_wtime()
  write(*,*) 'real time consumed for one iteration:',ompt2-ompt1,'seconds'
  write(*,*) 'cpu time consumed for one iteration:',t2-t1,'seconds'
enddo
if (flag_check) call checkpoint
write(*,*) 'END MULTIGRID'

write(*,*) ' '
write(*,*) 'def:'
write(*,*) 'avg=',real(sum(real(def,8))/float(ng1)/float(ng2)/float(ng3))
write(*,*) 'sigma=',sqrt(real(sum(real(def**2,8))/float(ng1)/float(ng2)/float(ng3)))
if (look) write(*,*) def(:,1,1)

!begin output deformation potential field
if (flag_def) then
  filename=trim(workdir)//trim(fndef)//trim(zstring)//'def'//trim(idstring)//gstring//trim(prefix)
  write(*,*) 'writing: ',trim(filename)
  open(31,file=filename,status='replace',access='stream')
  write(31) def
  close(31)
endif
!end output deformation potential field

!begin output gradient field
if (flag_gradient) then
  call cpu_time(t1)
  call scalar2gradient(def,grad,(/ng1,ng2,ng3/))

  filename=trim(workdir)//trim(fndisp)//trim(zstring)//'psix'//trim(idstring)//gstring//trim(prefix)
  write(*,*) 'writing: ',trim(filename)
  open(31,file=filename,status='replace',access='stream')
  write(31) grad(:,:,:,1)
  close(31)

  filename=trim(workdir)//trim(fndisp)//trim(zstring)//'psiy'//trim(idstring)//gstring//trim(prefix)
  write(*,*) 'writing: ',trim(filename)
  open(31,file=filename,status='replace',access='stream')
  write(31) grad(:,:,:,2)
  close(31)

  filename=trim(workdir)//trim(fndisp)//trim(zstring)//'psiz'//trim(idstring)//gstring//trim(prefix)
  write(*,*) 'writing: ',trim(filename)
  open(31,file=filename,status='replace',access='stream')
  write(31) grad(:,:,:,3)
  close(31)
  call cpu_time(t2)

  if (look) write(*,*) ' '
  write(*,*) 'time consumed for defp2phi:',t2-t1,'seconds'
  if (look) write(*,*) 'dx:'
  if (look) write(*,fmtstring) grad(:,1,1,1)
endif
! end output gradient field

! begin output reconstructed density field
if (flag_density) then
  call cpu_time(t1)
  call defp2delta(def,delta,(/ng1,ng2,ng3/))

  filename=trim(workdir)//trim(fndelta)//trim(zstring)//'recon'//trim(idstring)//gstring//trim(prefix)
  write(*,*) 'writing: ',trim(filename)
  open(31,file=filename,status='replace',access='stream')
  write(31) delta(:,:,:)
  close(31)

  write(*,*) ' '
  call cpu_time(t2)
  write(*,*) 'time consumed for defp2delta:',t2-t1,'seconds'
  write(*,*) 'delta:'
  if (look) write(*,fmtstring) delta(:,1,1)
  write(*,*) 'avg=',real(sum(real(delta(:,:,:),8))/float(ng1)/float(ng2)/float(ng3))
  write(*,*) 'sigma=',real(sqrt(sum(real(delta(:,:,:)**2,8))/float(ng1)/float(ng2)/float(ng3)))
endif
! end output reconstructed density field

contains

subroutine readrestart
implicit none
write(*,*) ' '
write(*,*) 'RESTART!'
write(sstring,'(i4.4)') restartstep
filename=trim(workdir)//trim(fnu)//trim(zstring)//trim(rstring_u)//sstring//'-'//trim(idstring)//gstring//trim(prefix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='old',access='stream')
read(31) u(1,:,:,:)
close(31)
filename=trim(workdir)//trim(fndef)//trim(zstring)//trim(rstring_d)//sstring//'-'//trim(idstring)//gstring//trim(prefix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='old',access='stream')
read(31) def
close(31)
write(*,*) 'def sigma=',real(sqrt(sum(real((def)**2,8))/float(ng1)/float(ng2)/float(ng3)))
write(*,*) 'u sigma=',real(sqrt(sum(real((u(1,:,:,:)-1.)**2,8))/float(ng1)/float(ng2)/float(ng3)))
endsubroutine readrestart

subroutine checkpoint
implicit none
write(*,*) ' '
write(*,*) 'CHECKPOINT!'
write(sstring,'(i4.4)') niter
filename=trim(workdir)//trim(fnu)//trim(zstring)//trim(rstring_u)//sstring//'-'//trim(idstring)//gstring//trim(prefix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) u(1,:,:,:)
close(31)
filename=trim(workdir)//trim(fndef)//trim(zstring)//trim(rstring_d)//sstring//'-'//trim(idstring)//gstring//trim(prefix)
write(*,*) 'reading: ',trim(filename)
open(31,file=filename,status='replace',access='stream')
write(31) def
close(31)
endsubroutine checkpoint

endprogram main
