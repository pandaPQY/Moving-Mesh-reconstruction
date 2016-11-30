! -*- Fortran -*-
! File Seidel.f
! written May 1994 by Ue-Li Pen, upen@astro.princeton.edu
! modified Nov 17 1995 to use static memory      
! We now use the global header files to allocate the biggest
! necessary amount of memory for the bottom level grid, and
! simply only use a smaller fraction of it at higher levels.
!    
#ifdef _T3D
#define BAJ6 8
#else
#define BAJ6 6
#endif  
#ifdef _SX5
! the extra padding reduces bank conflicts.
#define BAJI2 9+BAJ6
#else
#define BAJI2 8+BAJ6      
#endif
! for test runs uncomment next two lines
#ifdef TESTRELAX    
      call testrelax
      end
#endif
!
! Compiler options note:  on SGI, compile with
!      -WK,-ndr,-directives=-AK,-roundoff=2
!     so that it will ignore the convex and cray directives.
!     
!
! general styles and paradigms:
!
! temporary arrays:
!     are all allocated using the index of a common block stub.
! The interface wrapper then calls the actual work routine using
! the appropriate sectioning of the temporary arrays.
! This will result in an index violation at run time.  If compiling with
! array check, one needs to turn index checking off the the wrappers.
!
! Arguments:  I always try to put them in the order:
!  0.  everything < Temporary arrays      
!  1.  Arrays < Scalars
!  2.  intent(out) < intent(inout) < intent(in)
!  3.  floats < integers < logicals
!  4.  for arrays: bigger rank < smaller rank
!  5.  Arrays:     bigger size < smaller size
!  6.  For integers: dimensions < other
!  7.  Dimensions: same order as the arrays they declare
!
! Limitations:
!  the number of rows (na3) must be ab even multiple of nrelax.
!     
!  Machine specifics:
!    I have inserted vectorization and parallelization directives
!     for convex fc-8.0.  The algorithm is fully data parallel, so
!     it should run at full speed on a parallel/vector machine.
!    The convex compilers have some peculiar bugs, which parallelization
!     of the second layer loops due to some apparent varying inner trip
!     count.  This is fixed by vectorizing the second and parallelizing
!     the inner one.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!               
subroutine relax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu,nrelax,iopt)
implicit none
integer iopt,na1,na2,na3,nrelax,nu
real arr(na1,na2,na3), resid((na1+1)/2,(na2+1)/2,(na3+1)/2),&
     rhs(na1,na2,na3), rinf, u(nu,na1,na2,na3),&
     deformpot(na1,na2,na3),dx
!!dir$ shared *arr(:block,:block,:),*deformpot(:block,:block,:)
!!dir$ shared *rhs(:block,:block,:),*resid(:block,:block,:)
!!dir$ shared *u(:,:block,:block,:)
!
!
! perform a four color Gauss-Seidel relaxation given only the
! deformation potential at some reduced grid scale.  Optionally compute
! residual.
!
!
! the idea is to use a 2-D work array, so we maintain parallel/vector
! efficiency.
! Parameters:
!       arr(na1,na2,na3)                        the unknown field
!       deformpot(na1,na2,na3)   deformation potential
!       nrelax                           number of Gauss-Seidel iterations.
!                                        zero is legal.
!       resid(na1,na2,na3)              optional residual array
!       dx                              delta x, the grid spacing 
!       integer iopt                    the operation to be performed:
!   mod(iopt,8): 0 -- dont calculate residual
!                1 -- calculate residual on half grid size
!                2 -- return infinity norm of residual in rinf
!                3 -- return infinity norm in rinf, and residual
!                            in resid(,,)
!     The residual has the same sign as RHS, i.e. rhs-L[u].
!      
!   Let j=mod(iopt/8,8).  Then the procedure does one of the following for j:
!       1: regular Poisson iteration using metric g_{ab}
!       2: potential flow solver using triad e^a_b
!       3: implicit hydro solver (not implemented).
!
!  Return codes (passed in iopt):
!       iopt>=0:        successfully completed relaxation
!       iopt = -1:      invalid parameter
!
! it is assumed that all indeces integer powers of 2.
!
! ANSI Fortran-77 issues: I use ENDDO, long variable names, INCLUDE
! and IMPLICIT NONE.
! I have tried to make the variable names unique, even if truncated
! to six characters.  I do not know of compilers where it fails.
! The IMPLICIT NONE statement can simply be removed if not supported.
! The ENDDOs would need to be replaced by continue statements.
!
! parallel/porting issues: temporary arrays are allocated on a common
! block.  To be more portable, one could simply eliminate the common lines.
!
! up equal or less space than reals.
!      
! levels are counted increasing down the pyramid:
!  1                  X
!  2                  X   X
!  3                  X X X X
! etc.      
!

! locals:
integer i,nrmax,iiopt, nrmleft,itmax


! there are problems when nrelax is either too small, there will be
! problems with index wraparound.
!      
nrmax = min(4,max(1,na3/4))
! somehow, it doesnt like nrelax=2 for na3=8, so simply surpress it.      
nrmleft=mod(nrelax,nrmax)


!      write(*,*)'RELAX:def=',deformpot
!      write(*,*)'RELAX:rhs=',rhs
!      write(*,*)'RELAX:arr=',arr


if (nrelax .gt. nrmax) then
! dont calculate residuals until the last sweep         
   iiopt=iopt-mod(iopt,8)
   itmax=nrelax/nrmax
   do i=1,itmax
      if (i .eq. itmax .and. nrmleft .eq. 0) iiopt=iopt
      call r1lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu,nrmax,iiopt)
   enddo
   if ( nrmleft .gt. 0) then
      call r1lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu,nrmleft,iopt)
   endif
else
   call r1lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu,nrelax,iopt)
endif
!      write(*,*)'RELAX(exit):arr=',arr
return
end


!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      
#ifdef _ALPHA
subroutine r1lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu,nrelax,iopt)
#include "dimen.f90"
real flat(BAJ6*(BAJI2)*NG*NG+8192)
real flat2(BAJ6*(BAJI2)*NG*NG+8192)
! locals
integer iwordsize,ioff,ichk,ioff2
logical firsttime,firsthalf
save firsttime,ioff,ioff2,firsthalf
data firsttime /.true./
data firsthalf /.true./

!dec$ alias get_page_size, "getpagesize"
integer get_page_size,ipagesize
external get_page_size

if (firsttime) then
   firsttime=.false.
ipagesize=get_page_size()
iwordsize=%loc(flat(2))-%loc(flat(1))
ioff=(ipagesize-mod(%loc(flat(1)),ipagesize))/iwordsize+1
ioff2=(ipagesize-mod(%loc(flat2(1)),ipagesize))/iwordsize+1
write(*,*) 'r1lax: ioff=',ioff,ioff2
!dec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(flat)
endif
if (firsthalf .and. na1 .eq. NG/2) then
   firsthalf=.false.
   write(*,*) 'aligning halfsize baj',na1,na2,na3
!dec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(flat2)
endif

      
if (na1 .eq. NG) then
  call r2lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu,nrelax,iopt,flat(ioff))
else
  call r2lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu,nrelax,iopt,flat2(ioff2))
endif
return
end

subroutine r2lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu,nrelax,iopt,baj)
#else
subroutine r1lax(arr,resid,rhs,deformpot,u,rinf,dx,na1,na2,na3,nu,nrelax,iopt)
#endif
implicit none
integer iopt,na1,na2,na3,nrelax,iclvls,it,nu
real arr(na1,na2,na3), resid((na1+1)/2,(na2+1)/2,(na3+1)/2),&
     rhs(na1,na2,na3),rinf,u(nu,na1,na2,na3),&
     deformpot(na1,na2,na3),dx
!!dir$ shared *arr(:block,:block,:),*deformpot(:block,:block,:)
!!dir$ shared *rhs(:block,:block,:),*resid(:block,:block,:)
!!dir$ shared *u(:,:block,:block,:)
!
! Locals
!
#include "dimen.f90"
integer idxf(NG)
#define MAXRELAX 4
! I think I actually only use baj(1:2*nrelax+1)
      
#ifdef DYNAMIC_MEMORY
real baj(BAJ6,BAJI2,na1,na2)
!!dir$ shared baj(:,:,:block,:block)
#else      
!real baj(BAJ6,2*(MAXRELAX+1/(MAXRELAX+1))+BAJ6,NG,NG)
real baj(BAJ6,BAJI2,NG,NG)
#endif      
logical lsquare
integer k,nfudge,ki,ko,idateline,kp,iresidopt

#ifdef _ALPHA
!cdec$ MIGRATE_NEXT_TOUCH_NOPRESERVE(baj)
#endif

! always a special case for nrelax=0
nfudge=2*(nrelax+1/(nrelax+1))+6
nfudge=min(nfudge,na3)
if (iopt .lt. 0 .or. iopt .gt. 256 ) then
   write(*,*) 'relax: iopt exceeds legal range:',iopt
   pause
endif


if (mod(na3,2) .ne. 0) then
   write(*,*)'relax: dimensions must be even'
!   tmp=0
!   tmp=1/tmp
   pause
endif

if (mod(iopt/8,8) .eq. 1) then
   lsquare=.true.
else
   lsquare=.false.
endif
rinf=0
iresidopt=mod(iopt,8)
if (iresidopt .eq. 4) rinf=1.0e20
      
! computation strategy:
! the basic idea is red-black Gauss-Seidel.
!  There is a significant operation count to construct the triad or
!  metric, so we need to retain it for each of the sweeps.
!  Temporaries are kept in planes
!
! this looks rather complicated.  To illustrate an example
! with nrelax=2:
!
!    18 20  6 10  9 14 13 17 16 19
!     1  3  2  5  4  8  7 12 11 15
!----------------------------------
!     1  2  3  4  5  6  7  8  9 10     
!
! or with nrelax=3  (which we will use in the examples)
!
!      
!    44 47 46 48 15 21 20 27 26 33 32 38 37 42 41 45
!    39 43  6 10  9 14 13 19 18 25 24 31 30 36 35 40
!     1  3  2  5  4  8  7 12 11 17 16 23 22 29 28 34
!---------------------------------------------------
!     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
!
!     which unraveled is the sequence:
! 1 3 2 5 4 3 7 6 5 4 9 8 7 6 5 11 10 9 8 7 6 13 12 11 10 9 8 15 14 13 12 11 10
! 16 15 14 13 12 1 16 15 14 2 1 16 3 2 4
!      
if (iresidopt .eq. 1 .or. iresidopt .eq. 3) then
   call szero(resid,(na1+1)/2,(na2+1)/2,(na3+1)/2)
endif


! Now relax first planes

idxf(na3-1)=nfudge-1
idxf(na3)=nfudge
      
! first set up metric matrix for 3 tiers:
!  do the bottom plane.      

idateline=1
idxf(1)=idateline
idateline=mod(idateline,nfudge)+1
idxf(2)=idateline
! note that calcbajxy must be called with an odd plane
! plow na3,1      
call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge,1,lsquare)

if (nrelax .eq. 0) then
   if (iresidopt .eq. 0) then
      write(*,*) 'relax: nothing to do!'
   endif
   do k=2,na3,2
      idateline=mod(idateline,nfudge)+1
      kp=mod(k,na3)+1
      idxf(kp)=idateline
      idateline=mod(idateline,nfudge)+1
      kp=mod(kp,na3)+1
      idxf(kp)=idateline
      kp=mod(k,na3)+1
      call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge,kp,lsquare)
      call calcresid(resid,baj,arr,rhs,idxf,rinf,dx,na1,na2,na3,nfudge,k,iresidopt)
   enddo
   return
endif


idateline=mod(idateline,nfudge)+1
idxf(3)=idateline
idateline=mod(idateline,nfudge)+1
idxf(4)=idateline
! and 2:3      
call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge,3,lsquare)



! the initial pyramid:
!      
!                15
!           6 10  9 14 13
!     1  3  2  5  4  8  7 12 11
!---------------------------------------------------
!     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
!
! so this is #1: (1)
k=1
call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)

do ko=2,(nrelax-1)*4+1,2
   iclvls=(ko-2)/4
   idateline=mod(idateline,nfudge)+1
   kp=mod(ko+2,na3)+1
   idxf(kp)=idateline
   idateline=mod(idateline,nfudge)+1
   kp=mod(ko+3,na3)+1
   idxf(kp)=idateline
   kp=mod(ko+2,na3)+1
   call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge,kp,lsquare)
! to visualize: when ko=2, iclvls=0, and we do #2, #3 ( 3 2 )
! second time around, ko=4, iclvls=0, do       #4, #5 ( 5 4 )
! and the third time, ko=6, iclvls=1, do       #7 #8 #9 #10 (7 6 5 4)
! and so on.         
   do ki=0,iclvls
      k=mod(ko-2*ki+4*na3,na3)+1
      call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
      k=mod(ko-2*ki+4*na3-1,na3)+1
      call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
   enddo
! on the second time, fill #6 (3), etc 
   if (mod(ko,4) .eq. 0) then
      k=mod(ko/2,na3)+1
      call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
   endif
enddo

      
! do intermediate plane
!     
!                   21 20 27 26 33 32
!                         19 18 25 24 31 30
!                               17 16 23 22 29 28 
!---------------------------------------------------
!     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
!      
do ko=(nrelax-1)*4+2,na3-1,2
   k=ko
   kp=mod(ko+2,na3)+1
   idateline=mod(idateline,nfudge)+1
   idxf(kp)=idateline
   kp=mod(kp,na3)+1
   idateline=mod(idateline,nfudge)+1
   idxf(kp)=idateline
   k=mod(ko+2,na3)+1
   call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge,k,lsquare)
   do ki=1,nrelax*2,2
      k=ko-ki+2
      call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
      k=ko-ki+1
      call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
   enddo
   if (iresidopt .gt. 0 .and. ko .gt. (nrelax-1)*4+2 ) then
! we want tier 6 for nrelax=3
! and tier 2 for nrelax=1            
      k=mod(4*na3+ko-nrelax*2+1,na3)+1
      call calcresid(resid,baj,arr,rhs,idxf,rinf,dx,na1,na2,na3,nfudge,k,iresidopt)
   endif
enddo
      
! now clean up the leftovers: the Finale
!     
!    44 47 46 48                      38 37 42 41 45
!    39 43                                  36 35 40
!                                                 34
!---------------------------------------------------
!     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
!
!
if (.false.) then
idateline=mod(idateline,nfudge)+1
idxf(1)=idateline
idateline=mod(idateline,nfudge)+1
idxf(2)=idateline
! note that calcbajxy must be called with an odd plane
! plow na3,1      
call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge,1,lsquare)
endif
k=na3
! do #34 (16)
call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
! loop through #35 #36 #37 #38 (15 14 13 12)
do ki=2,nrelax
   k=mod(4*na3-2*ki+2,na3)+1
   call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
   k=mod(4*na3-2*ki+1,na3)+1
   call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
enddo
if (iresidopt .gt. 0) then
   k=mod(4*na3-nrelax*2+1,na3)+1
   call calcresid(resid,baj,arr,rhs,idxf,rinf,dx,na1,na2,na3,nfudge,k,iresidopt)
endif
k=3
idateline=mod(idateline,nfudge)+1
idxf(k)=idateline
k=4
idateline=mod(idateline,nfudge)+1
idxf(k)=idateline
k=3
call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge,k,lsquare)
do ko=1,nrelax-1
! the first time we enter the ki loop, it=1, ko=1, and we do a full
! step sweep: #39 #40 #41 #42 (1 16 15 14)
! the second time round, it=2, ko=1, and the first step only has one
! brick left.
   k=mod(4*na3+2*ko+2,na3)+1
   idateline=mod(idateline,nfudge)+1
   idxf(k)=idateline
   k=mod(4*na3+2*ko+3,na3)+1
   idateline=mod(idateline,nfudge)+1
   idxf(k)=idateline
   k=mod(4*na3+2*ko+2,na3)+1
   call calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nfudge,k,lsquare)
   do it=1,2
      do ki=1,nrelax-ko
! we update three new rows each iteration, so we need to alternate
! one and two double calcbajxys.               
         if ( it .ne. 2 .or. ki .ne. 1 ) then
            k=mod(4*na3+2*ko-2*ki+2*it-2,na3)+1
            call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
         endif
         k=mod(4*na3+2*ko-2*ki+2*it-3,na3)+1
         call relaxplane(arr,baj,rhs,idxf,dx,na1,na2,na3,nfudge,k)
      enddo
   enddo
! end it=1,2         
   if (iresidopt .gt. 0) then
! update the remaining residuals at the top and bottom:
! normalize by: if nrelax=2, ko=1, we want k=na3
      k=mod(4*na3+4*ko-2*nrelax-1,na3)+1
      call calcresid(resid,baj,arr,rhs,idxf,rinf,dx,na1,na2,na3,nfudge,k,iresidopt)
      k=mod(4*na3+4*ko-2*nrelax+1,na3)+1
      call calcresid(resid,baj,arr,rhs,idxf,rinf,dx,na1,na2,na3,nfudge,k,iresidopt)
   endif
enddo
! end ko=1,nrelax
if (iresidopt .gt. 0) then
    k=mod(k+1,na3)+1
    if (nrelax .eq. 1) then
      k=2
    endif
    call calcresid(resid,baj,arr,rhs,idxf,rinf,dx,na1,na2,na3,nfudge,k,iresidopt)
 endif
!
! do you really believe that worked?
return
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      
! compute the residual on a coarsened grid
!
! when called on an even grid, it calculates the residual at ki/2,
! and side effects resid(ki/2+1).  It uses accesses arrays at levels
! ki-2, ki-1, ki, ki+1      
!
subroutine calcresid(resid,fudge,phi,rhs,idxf,rinf,dx,na1,na2,na3,nfudge,ki,iopt)

! WORKER routine.  This is a very expensive routine.
!
! operation count: 46.375 A  and 26.875 M for nred=2
! typically, iter=4, this is called every eight relax iterations,
! and the effective cost should be 9.0312 FLOP/element to be
! added to the 67 for planerelax.      
! In practice, most machines spend much more time than this
! requires.  An exception is the cray, which does this very efficiently.
!      
implicit none
integer na1,na2,na3,nfudge,idxf(na3),ki,iopt
real resid((na1+1)/2,(na2+1)/2,(na3+1)/2),phi(na1,na2,na3),&
     fudge(BAJ6,BAJI2,na1,na2),rhs(na1,na2,na3),dx,rinf
!!dir$ shared *rhs(:block,:block,:),*resid(:block,:block,:)
!!dir$ shared *fudge(:,:,:block,:block),*phi(:block,:block,:)
! locals
integer i,j,km,ko,kp,k,jp,jm,ip,im,ib,ibp,ibm,jb,jbp,jbm,kb,kbp
#ifdef DYNAMIC_MEMORY
real rt(na1,na2,2)
!!dir$ shared rt(:block,:block,:)
#else      
real rt(NG,NG,2)
#endif
real b11p,b11m,b22p,b22m,b33p,b33m,tdiag,t1,t2,di2,tresid,rinfmax,rinfmin
logical lnorminf,lresid
#ifdef T3D
real rinf1
!!dir$ shared rinf1
#endif

if (mod(ki,2) .eq. 1) then
    write(*,*) 'calcresid called with odd k'
!          di2=0
!          di2=1/di2
    pause
endif


lresid=.false.
lnorminf=.false.
if (mod(iopt,2) .eq. 1) lresid=.true.
if (mod(iopt/2,2) .eq. 1) then
   lnorminf=.true.
endif
di2=1/dx**2
#if   defined(_SGI_SOURCE) || defined(_SX5)
! the SGI does not know how to promote tests!
rinfmax=0
rinfmin=1.e10
! 7.2 f77
!*$* assert do(serial)
!$omp parallel default(private) firstprivate(ki,na1,na2,na3,di2,lresid)
!$omp& shared(fudge,idxf,phi,rhs,rt) reduction(max:rinfmax)
!$omp& reduction(min:rinfmin)  
do ko=1,2
   k=ki-2+ko
   kp=mod(k,na3)+1
   km=mod(k+na3-2,na3)+1
! help the convex:            
!$dir prefer_vector
! help the SGI:         
!c*$* assert do(concurrent)         
!$doacross local(i,j,TRESID,t1,t2,tdiag,b33m,b33p,b22m,b22p,b11m,b11p
!$&         ,im,ip,j,jp,jm), reduction(rinfmax,rinfmin)
!$omp do
   do j=1,na2
      jp=mod(j,na2)+1
      jm=mod(j+na2-2,na2)+1
!$dir prefer_parallel_ext            
!*$* assert do(serial)
      do i=1,na1
         ip=mod(i,na1)+1
         im=mod(i+na1-2,na1)+1
! calculate a symmetric Laplacian
! \partial_a G^{ab} \partial_b \phi
! where G^{ab} is a 3x3 symmetric matrix given by fudge.
! We build a symmetric discretization by splitting the tensor
! sum into the diagonal and offdiagonal pieces, t1 and t2.
! The diagonal stencil is as follows:
!               
!                      o
!                     b22p 
!               o b11m o b11p o
!                     b22m
!                      o
!
! we get (\phi_{i+1}-\phi_i)*b_{11}^p-(\phi_i-\phi_{i-1})*b_{11}^m
! and so on.
!     
         b11p=(fudge(1,idxf(k),ip,j)+fudge(1,idxf(k),i,j))/2
         b11m=(fudge(1,idxf(k),im,j)+fudge(1,idxf(k),i,j))/2
         b22p=(fudge(4,idxf(k),i,jp)+fudge(4,idxf(k),i,j))/2
         b22m=(fudge(4,idxf(k),i,jm)+fudge(4,idxf(k),i,j))/2
         b33p=(fudge(6,idxf(kp),i,j)+fudge(6,idxf(k),i,j))/2
         b33m=(fudge(6,idxf(km),i,j)+fudge(6,idxf(k),i,j))/2
! 6 A,  6 M
         tdiag=-(b11p+b11m+b22p+b22m+b33p+b33m)
         t1=(b11p*phi(ip,j,k)+b11m*phi(im,j,k)+b22p*phi(i,jp,k) &
           +b22m*phi(i,jm,k)+b33p*phi(i,j,kp)+b33m*phi(i,j,km))
! 10 A, 6 M              
! t2 contains the diagonal terms
               t2=(  (phi(ip,jp,k)-phi(im,jp,k))*fudge(2,idxf(k),i,jp) &
                    -(phi(ip,jm,k)-phi(im,jm,k))*fudge(2,idxf(k),i,jm) &
                    +(phi(ip,jp,k)-phi(ip,jm,k))*fudge(2,idxf(k),ip,j) &
                    -(phi(im,jp,k)-phi(im,jm,k))*fudge(2,idxf(k),im,j) &
                                                                       &
                    +(phi(i,jp,kp)-phi(i,jm,kp))*fudge(5,idxf(kp),i,j) &
                    -(phi(i,jp,km)-phi(i,jm,km))*fudge(5,idxf(km),i,j) &
                    +(phi(i,jp,kp)-phi(i,jp,km))*fudge(5,idxf(k),i,jp) &
                    -(phi(i,jm,kp)-phi(i,jm,km))*fudge(5,idxf(k),i,jm) &
                                                                       &
                    +(phi(ip,j,kp)-phi(im,j,kp))*fudge(3,idxf(kp),i,j) &
                    -(phi(ip,j,km)-phi(im,j,km))*fudge(3,idxf(km),i,j) &
                    +(phi(ip,j,kp)-phi(ip,j,km))*fudge(3,idxf(k),ip,j) &
                    -(phi(im,j,kp)-phi(im,j,km))*fudge(3,idxf(k),im,j) &
                   )/4
! 3*8-1=23 A, 13 M
         tresid=rhs(i,j,k)-(t1+t2+tdiag*phi(i,j,k))*di2
! 3 A, 2 M
         rinfmax=max(rinfmax,abs(tresid))
         if (lresid) rt(i,j,ko)=tresid
         rinfmin=min(rinfmin,abs(tresid))
      enddo
   enddo
enddo
!$omp end parallel
if (lnorminf) rinf=max(rinf,rinfmax)
if (iopt .eq. 4) rinf=min(rinf,rinfmin)
#else
! zpj: I parallelized the following section by divide the condition
! zpj:  sentences if (lnormin) and if (iopt .eq.4 ) into two separate
! zpj:  section. (1)  if (lnormin) and (2) if (iopt .eq.4 ) because 
! zpj: (1) -> not (2) and (2) -> not (1)

if (lnorminf) then
! 7.2 f77
!$omp parallel default(private) shared(rt,fudge,phi,rhs) &
!$omp firstprivate(idxf,na1,na2,na3,di2,ki,lresid) reduction(max:rinf)
!*$* assert do(serial)
do ko=1,2
   k=ki-2+ko
   kp=mod(k,na3)+1
   km=mod(k+na3-2,na3)+1
! help the convex:         
!$dir prefer_vector
!!dir$ doshared(j,i) on rhs(i,j,1)      
! 7.2 f77
!*$* assert do(concurrent)
!$omp do
   do j=1,na2
      do i=1,na1
!$dir prefer_parallel_ext            
         jp=mod(j,na2)+1
         jm=mod(j+na2-2,na2)+1
         ip=mod(i,na1)+1
         im=mod(i+na1-2,na1)+1
! calculate a symmetric Laplacian
! \partial_a G^{ab} \partial_b \phi
! where G^{ab} is a 3x3 symmetric matrix given by fudge.
! We build a symmetric discretization by splitting the tensor
! sum into the diagonal and offdiagonal pieces, t1 and t2.
! The diagonal stencil is as follows:
!               
!                      o
!                     b22p 
!               o b11m o b11p o
!                     b22m
!                      o
!
! we get (\phi_{i+1}-\phi_i)*b_{11}^p-(\phi_i-\phi_{i-1})*b_{11}^m
! and so on.
!     
         b11p=(fudge(1,idxf(k),ip,j)+fudge(1,idxf(k),i,j))/2
         b11m=(fudge(1,idxf(k),im,j)+fudge(1,idxf(k),i,j))/2
         b22p=(fudge(4,idxf(k),i,jp)+fudge(4,idxf(k),i,j))/2
         b22m=(fudge(4,idxf(k),i,jm)+fudge(4,idxf(k),i,j))/2
         b33p=(fudge(6,idxf(kp),i,j)+fudge(6,idxf(k),i,j))/2
         b33m=(fudge(6,idxf(km),i,j)+fudge(6,idxf(k),i,j))/2
! 6 A,  6 M
         tdiag=-(b11p+b11m+b22p+b22m+b33p+b33m)
         t1=(b11p*phi(ip,j,k)+b11m*phi(im,j,k)+b22p*phi(i,jp,k) &
           +b22m*phi(i,jm,k)+b33p*phi(i,j,kp)+b33m*phi(i,j,km))
! 10 A, 6 M              
! t2 contains the diagonal terms
         t2=(  (phi(ip,jp,k)-phi(im,jp,k))*fudge(2,idxf(k),i,jp) &
              -(phi(ip,jm,k)-phi(im,jm,k))*fudge(2,idxf(k),i,jm) &
              +(phi(ip,jp,k)-phi(ip,jm,k))*fudge(2,idxf(k),ip,j) &
              -(phi(im,jp,k)-phi(im,jm,k))*fudge(2,idxf(k),im,j) &
              +(phi(i,jp,kp)-phi(i,jm,kp))*fudge(5,idxf(kp),i,j) &
              -(phi(i,jp,km)-phi(i,jm,km))*fudge(5,idxf(km),i,j) &
              +(phi(i,jp,kp)-phi(i,jp,km))*fudge(5,idxf(k),i,jp) &
              -(phi(i,jm,kp)-phi(i,jm,km))*fudge(5,idxf(k),i,jm) &
              +(phi(ip,j,kp)-phi(im,j,kp))*fudge(3,idxf(kp),i,j) &
              -(phi(ip,j,km)-phi(im,j,km))*fudge(3,idxf(km),i,j) &
              +(phi(ip,j,kp)-phi(ip,j,km))*fudge(3,idxf(k),ip,j) &
              -(phi(im,j,kp)-phi(im,j,km))*fudge(3,idxf(k),im,j) &
             )/4
! 3*8-1=23 A, 12 M
         tresid=rhs(i,j,k)-(t1+t2+tdiag*phi(i,j,k))*di2
! 3 A, 2 M
! zpj: I send these two if ( if (lnorminf) before max and if (iopt
!               .eq. 4) before min ahead of the parallel unit) 
! zpj: the original command is:
! zpj:         rinf=max(rinf,abs(tresid)) 
! zpj:         if (lresid) rt(i,j,ko)=tresid
! zpj:         rinf=max(rinf,abs(tresid)

         rinf=max(rinf,abs(tresid)) 
         if (lresid) rt(i,j,ko)=tresid
      enddo
   enddo
enddo
!$omp end parallel
endif
if (iopt .eq. 4) then
! 7.2 f77
!$omp parallel default(private) shared(rt,fudge,phi,rhs) &
!$omp firstprivate(idxf,na1,na2,na3,di2,ki,lresid) reduction(min:rinf)
!*$* assert do(serial)
do ko=1,2
   k=ki-2+ko
   kp=mod(k,na3)+1
   km=mod(k+na3-2,na3)+1
! help the convex:         
!$dir prefer_vector
!!dir$ doshared(j,i) on rhs(i,j,1)      
! 7.2 f77
!*$* assert do(concurrent)
!$omp do
   do j=1,na2
      do i=1,na1
!$dir prefer_parallel_ext            
         jp=mod(j,na2)+1
         jm=mod(j+na2-2,na2)+1
         ip=mod(i,na1)+1
         im=mod(i+na1-2,na1)+1
! calculate a symmetric Laplacian
! \partial_a G^{ab} \partial_b \phi
! where G^{ab} is a 3x3 symmetric matrix given by fudge.
! We build a symmetric discretization by splitting the tensor
! sum into the diagonal and offdiagonal pieces, t1 and t2.
! The diagonal stencil is as follows:
!               
!                      o
!                     b22p 
!               o b11m o b11p o
!                     b22m
!                      o
!
! we get (\phi_{i+1}-\phi_i)*b_{11}^p-(\phi_i-\phi_{i-1})*b_{11}^m
! and so on.
!     
         b11p=(fudge(1,idxf(k),ip,j)+fudge(1,idxf(k),i,j))/2
         b11m=(fudge(1,idxf(k),im,j)+fudge(1,idxf(k),i,j))/2
         b22p=(fudge(4,idxf(k),i,jp)+fudge(4,idxf(k),i,j))/2
         b22m=(fudge(4,idxf(k),i,jm)+fudge(4,idxf(k),i,j))/2
         b33p=(fudge(6,idxf(kp),i,j)+fudge(6,idxf(k),i,j))/2
         b33m=(fudge(6,idxf(km),i,j)+fudge(6,idxf(k),i,j))/2
! 6 A,  6 M
         tdiag=-(b11p+b11m+b22p+b22m+b33p+b33m)
         t1=(b11p*phi(ip,j,k)+b11m*phi(im,j,k)+b22p*phi(i,jp,k) &
           +b22m*phi(i,jm,k)+b33p*phi(i,j,kp)+b33m*phi(i,j,km))
! 10 A, 6 M              
! t2 contains the diagonal terms
         t2=(  (phi(ip,jp,k)-phi(im,jp,k))*fudge(2,idxf(k),i,jp) &
              -(phi(ip,jm,k)-phi(im,jm,k))*fudge(2,idxf(k),i,jm) &
              +(phi(ip,jp,k)-phi(ip,jm,k))*fudge(2,idxf(k),ip,j) &
              -(phi(im,jp,k)-phi(im,jm,k))*fudge(2,idxf(k),im,j) &
              +(phi(i,jp,kp)-phi(i,jm,kp))*fudge(5,idxf(kp),i,j) &
              -(phi(i,jp,km)-phi(i,jm,km))*fudge(5,idxf(km),i,j) &
              +(phi(i,jp,kp)-phi(i,jp,km))*fudge(5,idxf(k),i,jp) &
              -(phi(i,jm,kp)-phi(i,jm,km))*fudge(5,idxf(k),i,jm) &
              +(phi(ip,j,kp)-phi(im,j,kp))*fudge(3,idxf(kp),i,j) &
              -(phi(ip,j,km)-phi(im,j,km))*fudge(3,idxf(km),i,j) &
              +(phi(ip,j,kp)-phi(ip,j,km))*fudge(3,idxf(k),ip,j) &
              -(phi(im,j,kp)-phi(im,j,km))*fudge(3,idxf(k),im,j) &
             )/4
! 3*8-1=23 A, 12 M
         tresid=rhs(i,j,k)-(t1+t2+tdiag*phi(i,j,k))*di2
! 3 A, 2 M
         if (lresid) rt(i,j,ko)=tresid
         rinf=min(rinf,abs(tresid))
      enddo
   enddo
enddo
!$omp end parallel
endif
#  ifdef T3D
! this is the standard multiple processor reduction:
!!dir$ master
rinf1=rinf
!!dir$ end master
if (lnorminf) then
!!dir$ atomic update
   rinf1=max(rinf1,rinf)
!!dir$ barrier
!!dir$ suppress
   rinf=rinf1
else if (iopt .eq. 4) then
!!dir$ atomic update
   rinf1=min(rinf1,rinf)
!!dir$ barrier
!!dir$ suppress
   rinf=rinf1
endif
#  endif
#endif
! the next section uses (27A+7M)/2/nred**2 operations per element
if (lresid) then
!$omp parallel default(private) firstprivate(na1,na2,na3,ki) &
!$omp shared(resid,rt)
   k=ki/2
   kp=mod(k,na3/2)+1
   kb=1
   kbp=2
!$dir prefer_vector      
!!dir$ doshared(j,i) on rt(i,j,1)      
! 7.2 f77 
!*$* assert do(concurrent)
!$omp do
   do j=1,(na2+1)/2
!*$* assert do(serial)
      do i=1,(na1+1)/2
         ib=2*i-1
         ibp=mod(ib,na1)+1
         ibm=mod(ib+na1-2,na1)+1
         jb=2*j-1
         jbp=mod(jb,na2)+1
         jbm=mod(jb+na2-2,na2)+1
         resid(i,j,k)=resid(i,j,k)+rt(ib,jb,kb)/8 &
              +(rt(ibp,jb,kb)+rt(ibm,jb,kb)+rt(ib,jbp,kb) &
             +rt(ib,jbm,kb) +rt(ib,jb,kbp))/16 &
              +(rt(ibp,jbp,kb)+rt(ibp,jbm,kb)+rt(ibm,jbp,kb) &
              +rt(ibm,jbm,kb)+rt(ibp,jb,kbp)+rt(ibm,jb,kbp) &
              +rt(ib,jbp,kbp)+rt(ib,jbm,kbp))/32 &
              +(rt(ibp,jbp,kbp)+rt(ibp,jbm,kbp)+rt(ibm,jbp,kbp) &
              +rt(ibm,jbm,kbp))/64
! 18 A 4 M,  but we divide by eight for nred=2
! the kbm part was taken care of in the previous iteration.
         resid(i,j,kp)=resid(i,j,kp)+rt(ib,jb,kbp)/16 &
              +(rt(ibp,jb,kbp)+rt(ibm,jb,kbp)+rt(ib,jbp,kbp) &
              +rt(ib,jbm,kbp))/32 &
              +(rt(ibp,jbp,kbp)+rt(ibp,jbm,kbp)+rt(ibm,jbp,kbp) &
              +rt(ibm,jbm,kbp))/64
! 9 A 3 M, ditto         
      enddo
   enddo
!$omp end parallel
endif
return
end
!
!      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine relaxplane(phi,fudge,rhs,idxf,dx,na1,na2,na3,nfudge,k)
!
! WORKER routine.  This is the single most expensive routine, typically
! accounting for about 80% of the total execution time.
!
implicit none
integer na1,na2,na3,k,nfudge,idxf(na3)
real phi(na1,na2,na3),rhs(na1,na2,na3)
real fudge(BAJ6,BAJI2,na1,na2),dx
!!dir$ shared *rhs(:block,:block,:),*phi(:block,:block,:)
!!dir$ shared *fudge(:,:,:block,:block)
!
! Operation count: Floating point - 41 Additions, 26 Multiplicatins, 1 Division
!   per element (=> 67 FLOP),    plus various integer stuff
!      
! This should be compared to 6 Additions and 1 Multiplication
!      for a Euclidian Poisson iteration.
!      
! locals
real dx2,t1,t2,tdiag,b11p,b11m,b22p,b22m,b33p,b33m
integer i,j,io,jo,ip,im,jp,jm,kp,km
#ifdef _SX5
 integer ij
#endif
#ifdef _ALPHA
!dec$ alias addr_to_rad, "_OtsAddrToRad"
integer  addr_to_rad,irad,irad1,irad2
external addr_to_rad
integer omp_get_thread_num,cpuid,cpu_get_rad
external omp_get_thread_num,cpuid,cpu_get_rad
#endif

! executable section
!c$omp parallel default(private) 
!c$omp& shared(fudge,phi,rhs,dx,na1,na2,na3,idxf,k)
dx2=dx**2
kp=mod(k,na3)+1
km=mod(k+na3-2,na3)+1

! tell the cray preprocessor that there is no data dependency in this routine.
!fpp$ nodepchk r
! Integer division means that io,jo loop to 1 unless na=1.
!fpp$ skip
!$dir scalar      
do io=0,1-1/na1
!fpp$ skip
!$dir scalar         
   do jo=0,1-1/na2
! on the SGI Mipspro 7.2 compiler, the no recurrence directive
! no longer exists.  There the doacross will still force parallelization.
!fpp$ select concur
!$dir no_recurrence, force_parallel_ext
!!dir$ doshared(j,i) on phi(i,j,1)      
!c*$* assert no recurrence(phi)
! 7.2 f77
!*$* assert do(concurrent)
#ifdef _SX5
!cdir nodep
       do ij=0,na2*na1/4-1
          i=mod(ij,na1/2)*2+1+io
          j=ij*2/na1
          j=j*2+1+jo
#else  
!$omp parallel do &
!$omp private(i,j,ip,im,jp,jm,b11p,b11m,b22p,b22m,b33p,b33m,tdiag,t1,t2) &
!$omp shared(k,idxf,phi,fudge,rhs,dx2,na1,na2,io,jo)
      do j=1+jo,na2,2
#ifdef _ALPHAXXXXX
      if (k .eq. 1 .and. na2 .ge. 128) then
         if(mod(j,na2/32) .eq. 1) then
         irad1=addr_to_rad(rhs(1,j,1))
         irad2=addr_to_rad(fudge(1,1,1,j))
         irad=cpu_get_rad()
         if ((irad1 .ne. irad) .or. (irad2 .ne. irad)) then
                 write(*,*) j,irad1,irad2,irad,na2
         endif
         endif
         endif
#endif
!fpp$ select vector
!$dir no_recurrence, force_vector
!!dir$ ivdep
         do i=1+io,na1,2
#endif
         jp=mod(j,na2)+1
         jm=mod(j+na2-2,na2)+1
            ip=mod(i,na1)+1
            im=mod(i+na1-2,na1)+1
            b11p=(fudge(1,idxf(k),ip,j)+fudge(1,idxf(k),i,j))/2
            b11m=(fudge(1,idxf(k),im,j)+fudge(1,idxf(k),i,j))/2
            b22p=(fudge(4,idxf(k),i,jp)+fudge(4,idxf(k),i,j))/2
            b22m=(fudge(4,idxf(k),i,jm)+fudge(4,idxf(k),i,j))/2
            b33p=(fudge(6,idxf(kp),i,j)+fudge(6,idxf(k),i,j))/2
            b33m=(fudge(6,idxf(km),i,j)+fudge(6,idxf(k),i,j))/2
! 6 A,  6 M
            tdiag=-(b11p+b11m+b22p+b22m+b33p+b33m)
        t1=(b11p*phi(ip,j,k)+b11m*phi(im,j,k)+b22p*phi(i,jp,k) &
          +b22m*phi(i,jm,k)+b33p*phi(i,j,kp)+b33m*phi(i,j,km))
! 10 A, 6 M              
! t2 contains the diagonal terms
        t2=(   (phi(ip,jp,k)-phi(im,jp,k))*fudge(2,idxf(k),i,jp) &
              -(phi(ip,jm,k)-phi(im,jm,k))*fudge(2,idxf(k),i,jm) &
              +(phi(ip,jp,k)-phi(ip,jm,k))*fudge(2,idxf(k),ip,j) &
              -(phi(im,jp,k)-phi(im,jm,k))*fudge(2,idxf(k),im,j) &
              +(phi(i,jp,kp)-phi(i,jm,kp))*fudge(5,idxf(kp),i,j) &
              -(phi(i,jp,km)-phi(i,jm,km))*fudge(5,idxf(km),i,j) &
              +(phi(i,jp,kp)-phi(i,jp,km))*fudge(5,idxf(k),i,jp) &
              -(phi(i,jm,kp)-phi(i,jm,km))*fudge(5,idxf(k),i,jm) &
              +(phi(ip,j,kp)-phi(im,j,kp))*fudge(3,idxf(kp),i,j) &
              -(phi(ip,j,km)-phi(im,j,km))*fudge(3,idxf(km),i,j) &
              +(phi(ip,j,kp)-phi(ip,j,km))*fudge(3,idxf(k),ip,j) &
              -(phi(im,j,kp)-phi(im,j,km))*fudge(3,idxf(k),im,j) &
              )/4
! 3*8-1=23 A, 13 M
               phi(i,j,k)=(rhs(i,j,k)*dx2-t1-t2)/tdiag
! 2 A, 1 M, 1 D  
#ifndef _SX5                
         enddo
#endif
      enddo
      !$omp end parallel do
   enddo
enddo
!c$omp end parallel
return
end

!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
! calculate the metric at level k-1 and k, interpolating if necessary.
!



subroutine calcbajxy(baj,deformpot,u,idxf,dx,na1,na2,na3,nu,nbaj,k,lsquare)
implicit none
integer nu,na1,na2,na3,k,idxf(na3),nbaj
real u(nu,na1,na2,na3)
real  baj(BAJ6,BAJI2,na1,na2),deformpot(na1,na2,na3),dx
!!dir$ shared *deformpot(:block,:block,:),*u(:,:block,:block,:)
!!dir$ shared *baj(:,:,:block,:block)
logical lsquare
!      
! locals

integer km

if (mod(k,2) .ne. 1) then
   write(*,*)' calcbajxy must be called on odd level'
   pause
endif
km=mod(k+na3-2,na3)+1
call cbajraw(baj,deformpot,u,dx,na1,na2,na3,nu,nbaj,idxf(k),k,lsquare)
call cbajraw(baj,deformpot,u,dx,na1,na2,na3,nu,nbaj,idxf(km),km,lsquare)
return
end
      
      
subroutine cbajraw(baj,def,u,dx,na1,na2,na3,nu,nbaj,ik,k,lsquare)
! WORKER routine.  The second most expensive call in the whole program.
implicit none
logical lsquare
integer nu,na1,na2,na3,k,ik,nbaj
real baj(BAJ6,BAJI2,na1,na2),def(na1,na2,na3),dx,u(nu,na1,na2,na3)
!!dir$ shared *def(:block,:block,:),*u(:,:block,:block,:)
!!dir$ shared *baj(:,:,:block,:block)
! calculate the metric in the plane k
!
! Operation count:  28+12 A, 38+18 M, 1 D  =>  97 FLOP/element
!      
! locals
integer kp,km,j,jp,jm,i,ip,im
real phixx,phiyy,phizz,phixy,phixz,phiyz,a11,a12,a13,a22,a23,a33,&
     b11,b12,b13,b22,b23,b33,di2,det,urg
logical singular
!      real dtmin
!      save dtmin
!      data dtmin /1./
#ifdef _SX5
 integer ij
#endif
#ifdef GMETRIC
      include 'relaxgl.f90'
      include 'gmetric.f90'

if (ng1.eq.na1 .and. ng2.eq.na2 .and. ng3.eq.na3) then
   urg=1
!$omp parallel do default(private) firstprivate(na1,na2,na3,k,ik,isquare)
!$omp& shared(gbaj,baj,u)
!$omp do
   do j=1,na2
     do i=1,na1
         det=gbaj(7,i,j,k)
         b11=gbaj(1,i,j,k)
         b12=gbaj(2,i,j,k)
         b13=gbaj(3,i,j,k)
         b22=gbaj(4,i,j,k)
         b23=gbaj(5,i,j,k)
         b33=gbaj(6,i,j,k)
         if (lsquare) then
         baj(1,ik,i,j)=(b11**2+b12**2+b13**2)*det
         baj(2,ik,i,j)=(b11*b12+b12*b22+b13*b23)*det
         baj(3,ik,i,j)=(b11*b13+b12*b23+b13*b33)*det
         baj(4,ik,i,j)=(b12**2+b22**2+b23**2)*det
         baj(5,ik,i,j)=(b12*b13+b22*b23+b23*b33)*det
         baj(6,ik,i,j)=(b13**2+b23**2+b33**2)*det
         else
#ifndef NORHORG
          kp=mod(k,na3)+1
          km=mod(k+na3-2,na3)+1
          jp=mod(j,na2)+1
          jm=mod(j+na2-2,na2)+1
          ip=mod(i,na1)+1
          im=mod(i+na1-2,na1)+1
          urg=u(1,im,jm,km)+u(1,im,jm,k)+u(1,im,jm,kp) &
             +u(1,im,j,km)+u(1,im,j,k)+u(1,im,j,kp) &
             +u(1,im,jp,km)+u(1,im,jp,k)+u(1,im,jp,kp) &
             +u(1,i,jm,km)+u(1,i,jm,k)+u(1,i,jm,kp) &
             +u(1,i,j,km)+u(1,i,j,k)+u(1,i,j,kp) &
             +u(1,i,jp,km)+u(1,i,jp,k)+u(1,i,jp,kp) &
             +u(1,ip,jm,km)+u(1,ip,jm,k)+u(1,ip,jm,kp) &
             +u(1,ip,j,km)+u(1,ip,j,k)+u(1,ip,j,kp) &
             +u(1,ip,jp,km)+u(1,ip,jp,k)+u(1,ip,jp,kp)
          urg=urg/27.
#endif
         baj(1,ik,i,j)=b11*urg
         baj(2,ik,i,j)=b12*urg
         baj(3,ik,i,j)=b13*urg
         baj(4,ik,i,j)=b22*urg
         baj(5,ik,i,j)=b23*urg
         baj(6,ik,i,j)=b33*urg
         endif
      enddo
    enddo
!$omp end parallel
else
#endif

di2=1/dx**2
kp=mod(k,na3)+1
km=mod(k+na3-2,na3)+1
singular=.false.
! the 7.2 f77 gets confused if the urg is left dangling in a conditional
! assignment.
urg=1
!$dir prefer_vector
!!dir$ doshared(j,i) on def(i,j,1)      
! 7.2 f77
!*$* assert do(concurrent)
#ifdef _SX5
do ij=0,na1*na2-1
   i=mod(ij,na1)+1
   j=ij/na1+1
#else
!$omp parallel do default(private) &
!$omp shared(na1,na2,na3,dx,k,ik,lsquare,def,baj,u,di2,kp,km)
do j=1,na2
   do i=1,na1
#endif
      jp=mod(j,na2)+1
      jm=mod(j+na2-2,na2)+1
      ip=mod(i,na1)+1
      im=mod(i+na1-2,na1)+1
      phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))*di2
      phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))*di2
      phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))*di2
! 3*( 2 A, 2 M )            
      phixy=(def(ip,jp,k)-def(im,jp,k)-def(ip,jm,k)+def(im,jm,k))*di2/4
      phiyz=(def(i,jp,kp)-def(i,jp,km)-def(i,jm,kp)+def(i,jm,km))*di2/4
      phixz=(def(ip,j,kp)-def(im,j,kp)-def(ip,j,km)+def(im,j,km))*di2/4
! 3*( 3 A, 1 M )            
      a11=1+phixx
      a12=phixy
      a13=phixz
      a22=1+phiyy
      a23=phiyz
      a33=1+phizz
! det is \sqrt{g}            
      det=a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2-a12**2*a33
! 7 A, 11 M, 1 D
#ifndef _SGI_SOURCE
! the SGI compiler has problems            
!            singular = singular .or. det .le. 0
#endif            
! these bij are actually \sqrt{g} e^i_\alpha            
      b11=(a22*a33-a23**2)
      b12=(a13*a23-a12*a33)
      b13=(a12*a23-a13*a22)
      b22=(a11*a33-a13**2)
      b23=(a12*a13-a11*a23)
      b33=(a11*a22-a12**2)
! 6*( 1 A, 2 M )
#ifndef _SGI_SOURCE
      if (det .le. 1.e-5) then
! if a coarse mesh point gets tangled, skip this point.
         b11=1
         b12=0
         b13=0
         b22=1
         b23=0
         b33=1
         det=1
      endif
#endif
      if (lsquare) then
! we want \sqrt{g} g^{\alpha\beta}               
         baj(1,ik,i,j)=(b11**2+b12**2+b13**2)/det
         baj(2,ik,i,j)=(b11*b12+b12*b22+b13*b23)/det
         baj(3,ik,i,j)=(b11*b13+b12*b23+b13*b33)/det
         baj(4,ik,i,j)=(b12**2+b22**2+b23**2)/det
         baj(5,ik,i,j)=(b12*b13+b22*b23+b23*b33)/det
         baj(6,ik,i,j)=(b13**2+b23**2+b33**2)/det
! 6*( 2 A, 4 M ) plus 1 D
      else
! this branch is not taken very often.
#ifndef NORHORG
         urg=u(1,im,jm,km)+u(1,im,jm,k)+u(1,im,jm,kp) &
            +u(1,im,j,km)+u(1,im,j,k)+u(1,im,j,kp) &
            +u(1,im,jp,km)+u(1,im,jp,k)+u(1,im,jp,kp) &
            +u(1,i,jm,km)+u(1,i,jm,k)+u(1,i,jm,kp) &
            +u(1,i,j,km)+u(1,i,j,k)+u(1,i,j,kp) &
            +u(1,i,jp,km)+u(1,i,jp,k)+u(1,i,jp,kp) &
            +u(1,ip,jm,km)+u(1,ip,jm,k)+u(1,ip,jm,kp) &
            +u(1,ip,j,km)+u(1,ip,j,k)+u(1,ip,j,kp) &
            +u(1,ip,jp,km)+u(1,ip,jp,k)+u(1,ip,jp,kp)
         urg=urg/27.
#endif
!               urg=u(1,i,j,k)
!               urg=1
!
         baj(1,ik,i,j)=urg*b11/det
         baj(2,ik,i,j)=urg*b12/det
         baj(3,ik,i,j)=urg*b13/det
         baj(4,ik,i,j)=urg*b22/det
         baj(5,ik,i,j)=urg*b23/det
         baj(6,ik,i,j)=urg*b33/det
      endif
#ifndef _SX5
   enddo
#endif
enddo
!$omp end parallel do
#ifdef GMETRIC
endif
#endif
return
end



subroutine rgzerosum(arr,def,dx,na1,na2,na3)
! WORKER routine.  expensive.
! set the mean of arr() to zero, by subtracting a constant multiple
! of \sqrt{g}.
!
implicit none
integer na1,na2,na3,k
real arr(na1,na2,na3),dx ,def(na1,na2,na3)
!!dir$ shared *def(:block,:block,:), *arr(:block,:block,:)
! locals
integer kp,km,j,jp,jm,i,ip,im,ns1,ns2,ns3,ndim,ib,jb,kb, &
        io,ii,jo,ji,ko,ki,kbm,jbm,ibm
real phixx,phiyy,phizz,phixy,phixz,phiyz,a11,a12,a13,a22,a23,a33, &
     di2,det,dfact,cfact
real(8) dsum,asum,csum

#ifdef T3D
real asum1, dsum1
!!dir$ shared asum1,dsum1
intrinsic sum
#endif
#ifdef GMETRIC
include 'relaxgl.f90'
include 'gmetric.f90'

dsum=0.d0
if (ng1.eq.na1 .and. ng2.eq.na2 .and. ng3.eq.na3) then
   !$omp parallel shared(gbaj,na1,na2,na3) private(i,j,k) reduction(+:dsum)
   do k=1,na3
      !$omp do
      do j=1,na2
         do i=1,na1
            dsum=dsum+gbaj(7,i,j,k)
         enddo
      enddo
   enddo
   !$end parallel
else
#endif    

di2=1/dx**2
dsum=0.d0
#ifdef T3D
asum=sum(real(arr,8))
#else
csum=0
! 7.2 f77
!*$* assert do(serial)
!$omp parallel default(private) firstprivate(na1,na2,na3,di2) &
!$omp shared(def) reduction(+:dsum)      
do k=1,na3
   kp=mod(k,na3)+1
   km=mod(k+na3-2,na3)+1
!$dir prefer_vector      
!!dir$ doshared(j,i) on arr(i,j,1)      
! 7.2 f77
!*$* assert do(concurrent)
!$omp do
   do j=1,na2
      do i=1,na1
      jp=mod(j,na2)+1
      jm=mod(j+na2-2,na2)+1
         ip=mod(i,na1)+1
         im=mod(i+na1-2,na1)+1
         phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))*di2
         phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))*di2
         phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))*di2
         phixy=(def(ip,jp,k)-def(im,jp,k)-def(ip,jm,k)+def(im,jm,k))*di2/4
         phiyz=(def(i,jp,kp)-def(i,jp,km)-def(i,jm,kp)+def(i,jm,km))*di2/4
         phixz=(def(ip,j,kp)-def(im,j,kp)-def(ip,j,km)+def(im,j,km))*di2/4
         a11=1+phixx
         a12=phixy
         a13=phixz
         a22=1+phiyy
         a23=phiyz
         a33=1+phizz
         det=(a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2-a12**2*a33)
! This csum business only works if we start off with a regular
! initial grid.  The implicit assumption is that a distorted grid
! implies a non-linear region, where the exact form of the
! zero offset does not matter.  This is not true in the non-linear
! region.
!               if (det .gt. 1) csum=csum+det-1
! we use the same correction trick as in routine recut
         dsum=dsum+det
      enddo
   enddo
enddo
!$omp end parallel
#endif
#ifdef GMETRIC
endif
#endif
asum=0
di2=1/dx**2

#if 0
!!dir$ master
asum1=0
!!dir$ end master
!!dir$ atomic update      
asum1=asum1+asum
! now make sure  the sum reduction finished:
!!dir$ barrier
!!dir$ suppress
asum=asum1
#endif      
!      dsum=na1*na2*na3
!      if (abs(dsum-na1*na2*na3) .lt. csum) then
!         cfact=(dsum-na1*na2*na3)/csum
!      else
!         cfact=0
!      endif
      

!*$* assert do(serial)
!$omp parallel default(private) firstprivate(na1,na2,na3) &
!$omp shared(arr) reduction(+:asum)
do k=1,na3
! 7.2 f77
!*$* assert do(concurrent)
!$omp do
   do j=1,na2
      do i=1,na1
         asum=asum+arr(i,j,k)
      enddo
   enddo
enddo
!$omp end parallel
dfact=asum/(dsum)
csum=0
!      write(*,*)'rgzerosum: dfact=',dfact

#ifdef GMETRIC
if (ng1.eq.na1 .and. ng2.eq.na2 .and. ng3.eq.na3) then
!*$* assert do(serial)
!$omp parallel default(private) shared(gbaj,arr) &
!$omp reduction(+:csum) firstprivate(na1,na2,na3,dfact)
   do k=1,na3
! 7.2 f77
!*$* assert do(concurrent)
!$omp do
      do j=1,na2
         do i=1,na1
            det=gbaj(7,i,j,k)
            arr(i,j,k)=arr(i,j,k)-dfact*det
            csum=csum+det
         enddo
      enddo
   enddo
!$omp end parallel
else
#endif
!*$* assert do(serial)
!$omp parallel default(private) firstprivate(na1,na2,na3,di2,dfact) &
!$omp shared(arr,def) reduction(+:csum)
do k=1,na3
   kp=mod(k,na3)+1
   km=mod(k+na3-2,na3)+1
         
!$dir prefer_vector      
!!dir$ doshared(j,i) on def(i,j,1)      
! 7.2 f77
!*$* assert do(concurrent)
!$omp do
   do j=1,na2
      do i=1,na1
      jp=mod(j,na2)+1
      jm=mod(j+na2-2,na2)+1
         ip=mod(i,na1)+1
         im=mod(i+na1-2,na1)+1
         phixx=(def(ip,j,k)-2*def(i,j,k)+def(im,j,k))*di2
         phiyy=(def(i,jp,k)-2*def(i,j,k)+def(i,jm,k))*di2
         phizz=(def(i,j,kp)-2*def(i,j,k)+def(i,j,km))*di2
         phixy=(def(ip,jp,k)-def(im,jp,k)-def(ip,jm,k)+def(im,jm,k))*di2/4
         phiyz=(def(i,jp,kp)-def(i,jp,km)-def(i,jm,kp)+def(i,jm,km))*di2/4
         phixz=(def(ip,j,kp)-def(im,j,kp)-def(ip,j,km)+def(im,j,km))*di2/4
         a11=1+phixx
         a12=phixy
         a13=phixz
         a22=1+phiyy
         a23=phiyz
         a33=1+phizz
         det=(a11*a22*a33+2*a12*a13*a23-a13**2*a22-a11*a23**2-a12**2*a33)
!         if (det .gt. 1) det=det-cfact*(det-1)
         arr(i,j,k)=arr(i,j,k)-dfact*det
         csum=csum+det
      enddo
   enddo
enddo
!$omp end parallel
#ifdef GMETRIC
endif
#endif
if (abs(csum-na1*na2*na3) .gt. 0.5) then
!         write(*,*) 'rgzerosum: determinant fudge failed, csum=',csum
endif
return
end

      
subroutine zerosum(arr,n1,n2,n3)
! WORKER routine
implicit none
integer n1,n2,n3
real arr(n1,n2,n3)
!!dir$ shared  *arr(:block,:block,:)
! locals
integer i,j,k
real(8) sum1
#ifdef T3D
intrinsic sum

sum1=sum(real(arr,8))
#else
sum1=0.d0
!$omp parallel default(private) firstprivate(n1,n2,n3) shared(arr) &
!$omp reduction(+:sum1)
!*$* assert do(serial)
do k=1,n3
! 7.2 f77
!*$* assert do(concurrent)
!$omp do
   do j=1,n2
      do i=1,n1
         sum1=sum1+arr(i,j,k)
      enddo
   enddo
enddo
!$omp end parallel
#endif
sum1=sum1/(n1*n2*n3)
!*$* assert do(serial)
!$omp parallel default(private) shared(arr) firstprivate(n1,n2,n3,sum1)
do k=1,n3
!!dir$ doshared(j,i) on arr(i,j,1)      
! 7.2 f77
!*$* assert do(concurrent)
!$omp do
   do j=1,n2
      do i=1,n1
         arr(i,j,k)=arr(i,j,k)-sum1
      enddo
   enddo
enddo
!$omp end parallel
return
end


subroutine szero(a,nx,ny,nz)
! WORKER routine
implicit none
integer nx,ny,nz
real a(nx,ny,nz)
!!dir$ shared *a(:block,:block,:)
! locals
integer i,j,k


#ifdef T3D
a=0
#else      
!!dir$ doshared(j,i) on a(i,j,1)      
!$omp parallel private(i,j,k) firstprivate(nx,ny,nz) shared(a)
do i=1,nx
!$omp do
   do j=1,ny
      do k=1,nz
         a(i,j,k)=0
      enddo
   enddo
enddo
!$omp end parallel     
#endif      
return
end



      
!
!----------------------------------------------------------------------
! Now some self-test routines to test various aspects of the iteration.
!

subroutine testrelax
implicit none
!      call tst1drelax
call test3d
return
end

subroutine tst1drelax
implicit none
integer n1,n2,i,j,k,iopt,nrelax,it,ip,im
parameter(n1=16,n2=1)
real phi(n2,n2,n1), rhs(n2,n2,n1), dpot(n2,n2,n1), dx, resid(n1), rinf,exact(n1)
!$omp parallel default(private) shared(phi,rhs,dpot,exact)
do k=1,n1
!$omp do
   do j=1,n2
      do i=1,n2
         phi(i,j,k)=0
         rhs(i,j,k)=0
         dpot(i,j,k)=0
         exact(k)=0
      enddo
   enddo
enddo
!$omp end parallel

rhs(n2,n2,n1/2-1)=n1
call zerosum(rhs,n2,n2,n1)
iopt=0
dx=1
nrelax=2
do k=1,nrelax*2
do it=1,2
   do i=it,n1,2
      ip=mod(i,n1)+1
      im=mod(i+n1-2,n1)+1
      exact(i)=(exact(ip)+exact(im)+4*exact(i)-rhs(1,1,i))/6
   enddo
enddo
call zerosum(exact,n2,n2,n1)
enddo
 
call relax(phi,resid,rhs,dpot,rinf,dx,n2,n2,n1,nrelax,iopt)
call relax(phi,resid,rhs,dpot,rinf,dx,n2,n2,n1,nrelax,iopt)
call zerosum(phi,n2,n2,n1)
write(*,'(8f10.6)')(phi(1,1,j),j=1,n1)
write(*,*)'exact:'
write(*,'(8f10.6)')(exact(j),j=1,n1)
write(*,*)'real run:'
do i=1,100
   iopt=3
   call relax(phi,resid,rhs,dpot,rinf,dx,n2,n2,n1,nrelax,iopt)
   if (mod(i,10).eq. 0) then
      write(*,*) 'rinf=',rinf
      call zerosum(phi,n2,n2,n1)
      write(*,'(8f10.6)')(resid(j),j=1,n1)
      write(*,'(8f10.6)')(phi(1,1,j),j=1,n1)
   endif
enddo

iopt=2
nrelax=0
call relax(phi,resid,rhs,dpot,rinf,dx,n2,n2,n1,nrelax,iopt)
write(*,*) 'rinf=',rinf
return
end


subroutine test3d
implicit none
! Solve a 3-D periodic Poisson Problem, and compare with exact solution.
! Then measure the convergence rate, as function of nrelax, etc.
integer n
#include "dimen.f90"
parameter(n=NG)
real phi(n,n,n),rhs(n,n,n),dpot(n,n,n),exact(n,n,n),dx,rinf, resid(n,n,n), u(1,n,n,n)
integer n1,n2,i,j,k,iopt,nrelax

!$omp parallel default(private) shared(phi,rhs,dpot,resid,u)
do k=1,n
!$omp do
   do j=1,n
      do i=1,n
         phi(i,j,k)=0
         rhs(i,j,k)=0
         dpot(i,j,k)=0
         resid(i,j,k)=0
         u(1,i,j,k)=1
      enddo
   enddo
enddo
!$omp end parallel
rhs(1,1,1)=n**3
call zerosum(rhs,n,n,n)

iopt=3
nrelax=4
dx=1.
write(*,*) 'starting relax'
call relax(phi,resid,rhs,dpot,u,rinf,dx,n,n,n,1,nrelax,iopt)   
write(*,*) 'rinf=',rinf
      
return
end
