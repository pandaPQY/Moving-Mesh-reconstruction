! -*- Fortran -*-
!*********************************************************************
! 12 Oct 2016: extracted from drivers.fpp
!
!*********************************************************************


!****************************************************************
subroutine recut(def,tmp2,t1)
!      
! lay out a dense mesh over the specified region of the grid
!       
implicit none
include 'relaxgl.fi'
real tmp2(ng1+2,ng2,ng3),def(ng1,ng2,ng3),t1(ng1,ng2,ng3)

! locals
include 'globalpa.fi' 
real pi,wx,wy,wz,w000,w001,w010,&
     w011,w100,w101,w110,w111,sigma2, anorm, cosm_a, grow0, grow1,&
     xk,yk,zk,akernel,afreq,qk,pk,sumw, sum1,p1,&
     rho000,rho001,rho010,rho011,rho100,rho101,rho110,rho111
real phixx, phiyy, phizz, phixy, phixz, phiyz, a11, a12, a13, a22,&
     a23, a33, det, x, y, z, cfact, correctionvolume, detmin,&
     weight
integer icx,icy,icz,nbl,i,j,k,io,jo,ko, idist, ifold(ng3),&
     ii, iopt,i1,i2,i3,ip,im,jp,jm,kp,km, ix, iy, iz,&
     ixp, iyp, izp, iweight, nsample
parameter(nsample=3)
real(8)::dsum
real*4 densinit(ng1,ng2,ng3)
!mydist densinit(*,block,*)
common /densi/densinit

pi=4*atan(1.)
! detmin is the threshold above which we do any density corrections
detmin=.9
 
icrecut(1:3)=10  ! read from COSMOPA.DAT    
icx=icrecut(1)
icy=icrecut(2)
icz=icrecut(3)
write(*,91) icx,icy,icz
91   format(' recut: centering on ',3i5)

call definit(t1,rcut0,rcompress1,rcompress2)
! now move the regridded region onto the cluster:
dsum=0

!$omp parallel default(private) reduction(+:dsum)
!$omp& shared(def,t1,tmp2,icx,icy,icz,densinit) 
do k=1,ng3
!$omp do
   do j=1,ng2
      do i=1,ng1
! suffix "o" is source cell               
         io=i-ng1/2-icx
         jo=j-ng2/2-icy
         ko=k-ng3/2-icz
         io=mod(io+2*ng1-1,ng1)+1
         jo=mod(jo+2*ng2-1,ng2)+1
         ko=mod(ko+2*ng3-1,ng3)+1
         def(i,j,k)=t1(io,jo,ko)
         densinit(i,j,k)=0
         dsum=dsum+(tmp2(i,j,k)-1)**2
      enddo
   enddo
enddo
!$omp end parallel

write(*,*)'recut: intrinsic rms variation=',sqrt(real(dsum)/(ng1*ng2*ng3))
dsum=0
correctionvolume=0
!$omp parallel default(private) firstprivate(detmin)
!$omp& shared(def,tmp2,densinit,t1) reduction(+:sum,correctionvolume)
do k=1,ng3
!$omp do
   do j=1,ng2
      do i=1,ng1
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
         x=i+(def(ip,j,k)-def(im,j,k))/2
         y=j+(def(i,jp,k)-def(i,jm,k))/2
         z=k+(def(i,j,kp)-def(i,j,km))/2
         if (x .lt. 1) x=x+ng1
         if (y .lt. 1) y=y+ng2
         if (z .lt. 1) z=z+ng3
         ix=x
         iy=y
         iz=z
         wx=x-ix
         wy=y-iy
         wz=z-iz
         ix=mod(ix-1,ng1)+1
         iy=mod(iy-1,ng2)+1
         iz=mod(iz-1,ng3)+1
         w111=wx*wy*wz
         w110=wx*wy*(1-wz)
         w101=wx*(1-wy)*wz
         w100=wx*(1-wy)*(1-wz)
         w011=(1-wx)*wy*wz
         w010=(1-wx)*wy*(1-wz)
         w001=(1-wx)*(1-wy)*wz
         w000=(1-wx)*(1-wy)*(1-wz)
         ixp=mod(ix,ng1)+1
         iyp=mod(iy,ng2)+1
         izp=mod(iz,ng3)+1
         rho000=tmp2(ix,iy,iz)
         rho001=tmp2(ix,iy,izp)
         rho010=tmp2(ix,iyp,iz)
         rho011=tmp2(ix,iyp,izp)
         rho100=tmp2(ixp,iy,iz)
         rho101=tmp2(ixp,iy,izp)
         rho110=tmp2(ixp,iyp,iz)
         rho111=tmp2(ixp,iyp,izp)
         p1=rho000*w000+rho001*w001+rho010*w010&
           +rho011*w011+rho100*w100+rho101*w101&
           +rho110*w110+rho111*w111
         t1(i,j,k)=det 
         dsum=dsum+p1*det
! we want rg0sum to result in zero net fluctuations if p1=0.
         if (det .gt. detmin) correctionvolume=correctionvolume+(det-detmin)
         do ko=-nsample,nsample
            do jo=-nsample,nsample
               do io=-nsample,nsample
               iweight=(1+abs(ko/nsample))*(1+abs(jo/nsample))*(1+abs(io/nsample))
               weight=1./(2*nsample)**3/iweight
! trapezoidal rule integration over the original domain
               x=i+(def(ip,j,k)-def(im,j,k))/2+io/(2.*nsample)
               y=j+(def(i,jp,k)-def(i,jm,k))/2+jo/(2.*nsample)
               z=k+(def(i,j,kp)-def(i,j,km))/2+ko/(2.*nsample)
               if (x .lt. 1) x=x+ng1
               if (y .lt. 1) y=y+ng2
               if (z .lt. 1) z=z+ng3
               ix=x
               iy=y
               iz=z
               wx=x-ix
               wy=y-iy
               wz=z-iz
               ix=mod(ix-1,ng1)+1
               iy=mod(iy-1,ng2)+1
               iz=mod(iz-1,ng3)+1
               w111=wx*wy*wz
               w110=wx*wy*(1-wz)
               w101=wx*(1-wy)*wz
               w100=wx*(1-wy)*(1-wz)
               w011=(1-wx)*wy*wz
               w010=(1-wx)*wy*(1-wz)
               w001=(1-wx)*(1-wy)*wz
               w000=(1-wx)*(1-wy)*(1-wz)
               ixp=mod(ix,ng1)+1
               iyp=mod(iy,ng2)+1
               izp=mod(iz,ng3)+1
               rho000=tmp2(ix,iy,iz)
               rho001=tmp2(ix,iy,izp)
               rho010=tmp2(ix,iyp,iz)
               rho011=tmp2(ix,iyp,izp)
               rho100=tmp2(ixp,iy,iz)
               rho101=tmp2(ixp,iy,izp)
               rho110=tmp2(ixp,iyp,iz)
               rho111=tmp2(ixp,iyp,izp)
               p1=rho000*w000+rho001*w001+rho010*w010&
                 +rho011*w011+rho100*w100+rho101*w101&
                 +rho110*w110+rho111*w111
! we subtract one to be less prone to rounding error.
               densinit(i,j,k)=densinit(i,j,k)+weight*(p1*det-1)
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end parallel
cfact=(dsum-ng1*ng2*ng3)/correctionvolume
dsum=dsum/(ng1*ng2*ng3)
write(*,*) 'mean subgridded rho=',real(dsum),'  fudge factor=',cfact
! now we use the underdense region to correct the mass excess/deficit
dsum=0

!$omp parallel default(private) firstprivate(detmin,cfact)
!$omp& shared(t1,densinit,tmp2) reduction(+:sum)
do k=1,ng3
!$omp do
   do j=1,ng2
      do i=1,ng1
         det=t1(i,j,k)
         if ( det .gt. detmin) then
            densinit(i,j,k) = densinit(i,j,k)-cfact*(det-detmin)
         endif
! this trading makes us a bit less prone to rounding error.  tmp
! is usually double precision, while densinit is single precision.
         tmp2(i,j,k)=densinit(i,j,k)+1
         densinit(i,j,k)=t1(i,j,k)
         dsum=dsum+tmp2(i,j,k)
      enddo
   enddo
enddo
!$omp end parallel

dsum=dsum/(ng1*ng2*ng3)
write(*,*) 'corrected subgridded rho=',real(dsum)
!
return
endsubroutine recut
      
