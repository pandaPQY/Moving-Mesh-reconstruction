program from1024to128
   implicit none
   integer,parameter :: ng=1024
   integer,parameter :: n_coarse=3
   integer,parameter :: n_final=ng/2**n_coarse
   integer nf
   real den(ng,ng,ng)
   real den_final(n_final,n_final,n_final)
   integer i,j,k,n
   character(200),parameter :: inputfile='/scratch2/p/pen/qiaoyin/simulation_1024/node0/0.000den0_2.bin'
   character(200),parameter :: outputfile='0.000den128_2.bin'
   open(11,file=trim(adjustl(inputfile)),access='stream')
   print*,'reading ',trim(adjustl(inputfile))
   read(11) den
   close(11)
      nf=ng
   do i=1,n_coarse
      nf=nf/2
      call coarsen(den(1:2*nf,1:2*nf,1:2*nf),nf)
      print*,'nf = ',nf
!      den_final(1:nf,1:nf,1:nf)=den(1:nf,1:nf,1:nf)
   enddo
   den_final=den(1:nf,1:nf,1:nf)
   print*,'sum of den_final = ',sum(den_final),' mean of den_final = ',sum(den_final)/n_final**3
   print*,'maximum = ',maxval(den_final),' minimum = ',minval(den_final)
!   den_final=real(den_final,8)/sum(real(den_final,8))*float(nf)**3
!   print*,'sum of den_final = ',sum(den_final),' mean of den_final = ',sum(den_final)/n_final**3
!   print*,'maximum = ',maxval(den_final),' minimum = ',minval(den_final)
!   n=0
!   do k=1,nf
!   do j=1,nf
!   do i=1,nf
!      if ( den_final(i,j,k)<10e-1 ) then
!          den_final(i,j,k)=0.1
!          print*,i,j,k
!          n=n+1
!      endif
!   enddo
!   enddo
!   enddo
!   print*,'n = ',n
!   den_final=real(den_final,8)/sum(den_final)**n_final**3-1
   print*,'sum of den_final = ',sum(den_final),' mean of den_final = ',sum(den_final)/n_final**3
   print*,'maximum = ',maxval(den_final),' minimum = ',minval(den_final)

   open(11,file=trim(adjustl(outputfile)),access='stream')
   print*,'writing ',trim(adjustl(outputfile))
   write(11) den_final
   close(11)
end     

subroutine coarsen(den,nf)
  implicit none
  integer nf
  integer i,j,k
  real den(2*nf,2*nf,2*nf)
  do k=1,nf
  do j=1,nf
  do i=1,nf
     den(i,j,k)=sum(den(2*i-1:2*i,2*j-1:2*j,2*k-1:2*k))
  enddo
  enddo
  enddo
end subroutine coarsen
