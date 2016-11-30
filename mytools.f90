subroutine setzero(arr,ng1,ng2,ng3)
implicit none
integer::ng1,ng2,ng3
real::arr(ng1,ng2,ng3)
integer::k
!$omp parallel do default(shared) private(k)
do k=1,ng3
  arr(:,:,k)=0.
enddo
!$omp end parallel do
endsubroutine setzero

subroutine calcsigma(arr,mean,sigma,ng1,ng2,ng3,vname)
implicit none
integer::ng1,ng2,ng3
real::arr(ng1,ng2,ng3),mean,sigma,ssum
character(*)::vname
integer::k
ssum=0.
!$omp parallel do private(k) reduction(+:ssum)
do k=1,ng3
  ssum=sum(real((arr(:,:,k)-mean)**2,8))
enddo
!$omp end parallel do
sigma=sqrt(ssum/ng1/ng2/ng3)
write(*,*) trim(vname),' sigma=',sigma
endsubroutine calcsigma

subroutine calcmean(arr,mean,ng1,ng2,ng3,vname)
implicit none
integer::ng1,ng2,ng3
real::arr(ng1,ng2,ng3),mean,ssum
character(*)::vname
integer::k
ssum=0.
!$omp parallel do private(k) reduction(+:ssum)
do k=1,ng3
  ssum=sum(real(arr(:,:,k),8))
enddo
!$omp end parallel do
mean=ssum/ng1/ng2/ng3
write(*,*) trim(vname),' mean=',mean
endsubroutine calcmean

subroutine matadd(arr3,arr1,arr2,ng1,ng2,ng3)
implicit none
integer::ng1,ng2,ng3
real::arr3(ng1,ng2,ng3),arr1(ng1,ng2,ng3),arr2(ng1,ng2,ng3)
integer::k
!$omp parallel do
do k=1,ng3
  arr3(:,:,k)=arr1(:,:,k)+arr2(:,:,k)
enddo
!$omp end parallel do
endsubroutine matadd
