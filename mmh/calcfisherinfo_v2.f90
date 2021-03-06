program calcfishinfo
   implicit none
   integer,parameter :: n_sim=136
   integer,parameter :: n_version=4
   integer,parameter :: nk=31
   integer,parameter :: nkraw=512
   real,parameter :: pi=3.1415926536
   real power_all(n_version,n_sim,10,nk),tmp_pow(nk)
   real kbin(nk),error(n_version,nk),tmp_invcov(nk,nk)
   real pow(9,n_sim,nk),pow_fine(n_version,n_sim,nk),correlation(n_version,n_sim,nk),pow1(n_sim),pow2(n_sim),meanpow(n_version+1,nk),meancorrelation(n_version,nk)
   real cov(n_version,nk,nk),covnorm(n_version,nk,nk),fishinfo(n_version+1,nk),invcov(nk,nk),checkinv(nk,nk),corr(n_version,nk,nk),tmp_fishinfo(1,nk),tmpp_fishinfo(nk,1)
   integer i_sim,xbin,ybin,i,j,k,tmp_ps
   character(len=3) whichbin,whichsim
   fishinfo=0
   invcov=0
   cov=0
   covnorm=0
   pow1=0
   pow2=0
   error=0
   meanpow=0
   print*,'nk=',nk

!!! read in power spectrum of nonlinear density field, reconstructed density
!!! field, linear density field @ z =0 and @ z =100 and gaussian random field
   open(22,file='power_all_256_0.dat',access='stream')
   read(22) power_all
   close(22)
   meanpow(1,:)=power_all(1,1,2,:)! k in physical unit
!   print*,'check k mode: ',power_all(1,:,2,1)
   do k=1,n_version
      tmp_pow=0
      do i_sim=1,n_sim
         tmp_pow=tmp_pow+reshape(power_all(k,i_sim,3,:),(/nk/))!/power_all(5,i_sim,3,:)
      enddo
      meanpow(k+1,:)=tmp_pow/n_sim
   enddo
   meancorrelation=sum(power_all(:,:,8,:),2)/n_sim
   do xbin=1,nk
         fishinfo(n_version+1,xbin)=sum(power_all(1,1,1,:xbin))!+power_all(1,1,1,xbin)
   enddo
   print*,'number of k=',fishinfo(n_version+1,:)
   do k=1,n_version
      do xbin=1,nk
         call calerror(power_all(k,:,3,xbin),error(k,xbin),n_sim)
      enddo
   enddo

   print*,'read powerspectrum done.'
   print*,'errorbar = ',error
   print*,'meanpow=',meanpow
   print*,'meancorrelation=',meancorrelation
!   meancorrelation=1!
   do k=1,n_version
      do ybin=1,nk
      do xbin=1,nk
            pow1=reshape(power_all(k,:,3,xbin),(/n_sim/))
            pow2=reshape(power_all(k,:,3,ybin),(/n_sim/))
            call calccov(pow1,pow2,n_sim,cov(k,xbin,ybin),covnorm(k,xbin,ybin))
            if (xbin==ybin ) then
!               if (k==5) print*,'cov=',cov(:,xbin,xbin)
               call inverse(reshape(covnorm(k,1:xbin,1:ybin),(/xbin,ybin/)),invcov(1:xbin,1:ybin),xbin)
!               if (k==5) print*,'inv=',invcov(xbin,xbin)
          checkinv(1:xbin,1:ybin)=MATMUL(reshape(covnorm(k,1:xbin,1:ybin),(/xbin,ybin/)),invcov(1:xbin,1:ybin))
!          do i=1,xbin
!          print*,xbin,'check inv',covnorm(k,i,i),invcov(i,i),checkinv(i,i)
!          enddo
!               tmp_fishinfo=0
!               tmpp_fishinfo=0
!!              invcov(i,j)=invcov(i,j)*meancorrelation(k,i)*meancorrelation(k,j)!*meanpow(6,i)*meanpow(6,j)/meanpow(k+1,i)/meanpow(k+1,j)
!               do i=1,xbin
!!               do j=1,ybin
!                  fishinfo(k,xbin)=fishinfo(k,xbin)+invcov(i,i)*meancorrelation(k,i)
!               enddo
!               enddo
               tmp_invcov(1:xbin,1:ybin)=TRANSPOSE(invcov(1:xbin,1:ybin))
               tmp_fishinfo(1:1,1:ybin)=MATMUL(meancorrelation(k:k,1:xbin)**2,invcov(1:xbin,1:ybin))
               tmpp_fishinfo(1:xbin,1:1)=TRANSPOSE(meancorrelation(k:k,1:xbin)**2)
              
               fishinfo(k:k,xbin:xbin)=MATMUL(tmp_fishinfo(1:1,1:xbin),tmpp_fishinfo(1:xbin,1:1))
                                 !           tmp_fishinfo(1,xbin)=1./covnorm(k,xbin,xbin)
                                 !           fishinfo(k,xbin)=sum(tmp_fishinfo(1,1:xbin))
           endif
      enddo
      enddo
   enddo

! do k=1,5
!   xbin=nk
!      call inverse(reshape(covnorm(k,1:xbin,1:xbin),(/xbin,xbin/)),invcov(1:xbin,1:xbin),xbin)
!    do xbin=1,nk
!      do i=1,xbin
!      do j=1,xbin
!         fishinfo(k,xbin)=fishinfo(k,xbin)+invcov(i,j)!*meancorrelation(k,i)
!      enddo
!      enddo
!   enddo
!enddo


   do xbin=1,nk
      write(whichbin,'(i3)') xbin
      print*,'nb,rec,lin fisherinfo(k=',trim((adjustl(whichbin))),')=',fishinfo(:,xbin)
   enddo
   checkinv=MATMUL(reshape(covnorm(n_version,:,:),(/nk,nk/)),invcov)
   do i=1,nk
   print*,'check inv',invcov(i,i),checkinv(i,i)
   enddo
   do ybin=1,nk
   do xbin=1,nk
      do k=1,n_version
         corr(k,xbin,ybin)=cov(k,xbin,ybin)/sqrt(cov(k,xbin,xbin)*cov(k,ybin,ybin))
      enddo
   enddo
   enddo
!   print*,'corr = ',corr(1,:,:)
   print*, 'fishinfo written!'
   do i=1,nk
!      print*,invcov(i,i),checkinv(i,i)
   enddo
   open(24,file='meanpower_all_v256.dat',status='replace',access='stream')!v2: delete noice                 v3: gaussian random field times transfer function to be linear field
   open(25,file='fishinfo_all_v256.dat',status='replace',access='stream')! v2: set meancorrelation to be 1  v4: set meancorrelation to be 1
   open(26,file='errorbar_all_v256.dat',status='replace',access='stream')                                  !v6: diagonal only
   open(27,file='corr_all_v256.dat',status='replace',access='stream')                                      !v5: covariance=../n_sim
   open(28,file='meancorrelation_all_v256.dat',status='replace',access='stream')                           !v256: rec--> lin
   write(24) meanpow
   write(25) fishinfo
   write(26) error
   write(27) corr
   write(28) meancorrelation
   close(24)
   close(25)
   close(26)
   close(27)
   close(28)
end

subroutine calerror(pow,error,n_sim)
   implicit none
   integer n_sim,i
   real pow(n_sim),pow_tmp,error
   pow_tmp=0
   error=sum(pow)/n_sim
do i=1,n_sim
   pow_tmp=pow_tmp+(pow(i)-error)**2
enddo
   error=sqrt(pow_tmp/(n_sim-1))!/sqrt(real(n_sim))
endsubroutine calerror

subroutine calccov(pow1,pow2,n_sim,cov12,covnorm12)

   implicit none
   integer n_sim,i
   real pow1(n_sim),pow2(n_sim),cov12,covnorm12,pow11(n_sim),pow22(n_sim)
   pow11=pow1-sum(pow1)/n_sim
   pow22=pow2-sum(pow2)/n_sim
!   do i=1,n_sim
!      cov12=cov12+pow11(i)*pow22(i)
!   enddo
   cov12=DOT_PRODUCT(pow11,pow22)/(n_sim-1)
!   cov12=DOT_PRODUCT(pow1,pow2)/n_sim-sum(pow1)*sum(pow2)/n_sim**2
   covnorm12=cov12/sum(pow1)/sum(pow2)*n_sim**2
!   cov12=abs(cov12)
endsubroutine calccov

subroutine inverse(A,IA,N)
  implicit none
  integer N
  real A(N,N),IA(N,N)
  real,allocatable :: B(:,:)
  integer :: i,j
  allocate(B(N,N))
  forall(i=1:N,j=1:N,i==j) IA(i,j)=1.0
  forall(i=1:N,j=1:N,i/=j) IA(i,j)=0.0
  ! Save original matrix A and array B
  B=A
  ! Set B into a diagonal matrix
  call Upper(B,IA,N) ! Set B into an upper triangular matrix
  call Lower(B,IA,N) ! Set B into a lower triangular matrix
  ! Solve
  forall(i=1:N) IA(i,:)=IA(i,:)/B(i,i)
  return
end subroutine inverse

!! Output submatrix
!subroutine output(matrix)
!  implicit none
!  real  :: matrix(:,:)
!  integer :: m,n,i
!  character(len=20) :: FOR='(??(1x,f6.3))'
!  m = size(matrix,1)
!  n = size(matrix,2)
!  ! Set output format
!  write( FOR(2:3), '(I2)' ) N
!  do i=1,N
!        write( *, FMT=FOR ) matrix(i,:)
!  end do
!  return
!end subroutine output
!
!
! Give submatrix of the upper triangular matrix 
subroutine Upper(M,S,N)
  implicit none
  integer :: N
  real    :: M(N,N)
  real    :: S(N,N)
  integer :: I,J
  real :: E
  do I=1,N-1
    do J=I+1,N
      E=M(J,I)/M(I,I)
      M(J,I:N)=M(J,I:N)-M(I,I:N)*E
      S(J,:)=S(J,:)-S(I,:)*E
    end do
  end do
  return
end subroutine Upper
! Give submatrix of the lower triangular matrix 
subroutine Lower(M,S,N)
  implicit none
  integer :: N
  real    :: M(N,N)
  real    :: S(N,N)
  integer :: I,J
  real :: E
  do I=N,2,-1
    do J=I-1,1,-1
      E=M(J,I)/M(I,I)
      M(J,1:N)=M(J,1:N)-M(I,1:N)*E
      S(J,:)=S(J,:)-S(I,:)*E
    end do
  end do
  return
end subroutine Lower

