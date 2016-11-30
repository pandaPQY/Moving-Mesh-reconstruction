#include "dimen.f90"
! -*- Fortran -*-

integer ng1,ng2,ng3, maxfluidcmp
integer nstore, ncstore
parameter(nstore=5,ncstore=15,maxfluidcmp=5,ng1=NG,ng2=ng1,ng3=ng1)
real wslopelim, wbee
parameter(wslopelim=1.5,wbee=1)
! these two parameters are to choose a TVD limiter.  If either wbee=0
! or wslopelim=1, the scheme becomes the most conservative second order
! TVD scheme: minmod.  The other extreme of the stability limit is
! wslopelim=2, wbee=1, which is called the "superbee" limiter.
! The minmod limiter tends to be quite diffusive.
! Unfortunately, the relaxing TVD scheme is not very TVD by the time
! one reaches the superbee.  In 1-D pancake tests, one sees features
! in the solution near the onset of the compression limiter.
! A reasonable compromise appears to be wslopelim=1.5, wbee=1
