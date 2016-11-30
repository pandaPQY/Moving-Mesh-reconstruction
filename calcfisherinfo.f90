program calcfishinfo
   implicit none
   integer,parameter :: n_sim=136
   integer,parameter :: nbin=48 
   integer,parameter :: nk=nbin-6
   integer,parameter :: nkraw=512
   double precision,parameter :: pi=3.1415926536
   double precision tmp(9,n_sim,2,nbin),tmp_kbin(nbin),tmp_lin(n_sim,5,nkraw),tmp_lin100(n_sim,5,nkraw),tmp_fishinfo(nk)
   double precision mp_kbin(nbin),kbin(nk),error(5,nk)
   double precision pow(9,n_sim,nk),pow_fine(5,n_sim,nk),correlation(5,n_sim,nk),pow1(n_sim),pow2(n_sim),meanpow(5,nk),meanpow_fine(5,nk),meancorrelation(5,nk)
   double precision cov(5,nk,nk),covnorm(5,nk,nk),fishinfo(5,nk),invcov(nk,nk),checkinv(nk,nk),corr(5,nk,nk),stn(5,nk)
   integer i_sim,xbin,ybin,i,j,k,tmp_ps
   character(len=3) whichbin,whichsim
   fishinfo=0
   invcov=0
   cov=0
   covnorm=0
   pow1=0
   pow2=0
   error=0
   stn=0
   print*,'nk=',nk
   do i_sim=2,n_sim+1
      i=i_sim-1
      write(whichsim,'(i3)') i_sim
      open(22,file='ps48bins_128_new/ps48bins_128_'//trim((adjustl(whichsim)))//'.txt',form='formatted',action='read')
      open(23,file='ps48bins_def_128_new/ps48bins_def_128_'//trim((adjustl(whichsim)))//'.txt',form='formatted',action='read')
      open(24,file='ps48bins_lin0_128/ps48bins_lin0_'//trim((adjustl(whichsim)))//'.txt',form='formatted',action='read')
      open(25,file='ps48bins_lin100_128/ps48bins_lin100_'//trim((adjustl(whichsim)))//'.txt',form='formatted',action='read')
      open(26,file='ps_gaussianrdf/ps_gaussianrdf_'//trim((adjustl(whichsim)))//'.txt',form='formatted',action='read')
!      open(24,file='0.000pk.init_300Mpc/pk_'//trim((adjustl(whichsim)))//'.init',form='formatted',action='read')
!      open(25,file='100.000pk.init_300Mpc/pk_'//trim((adjustl(whichsim)))//'.init',form='formatted',action='read')
      open(27,file='crossps_linnonlin/crossps_linnonlin_'//trim((adjustl(whichsim)))//'.txt',form='formatted',action='read')
      open(28,file='crossps_lindef/crossps_lindef_'//trim((adjustl(whichsim)))//'.txt',form='formatted',action='read')
      open(29,file='crossps_linlin0/crossps_linlin0_'//trim((adjustl(whichsim)))//'.txt',form='formatted',action='read')
      open(30,file='crossps_linlin100/crossps_linlin100_'//trim((adjustl(whichsim)))//'.txt',form='formatted',action='read')
      read(22,*) tmp(1,i,:,:)!n-body
      read(23,*) tmp(2,i,:,:)!reconstuction
      read(24,*) tmp(3,i,:,:)!linear at z=0
!      read(24,*) tmp_lin(i,:,:)
      read(25,*) tmp(4,i,:,:)!linear at z=100
      read(26,*) tmp(5,i,:,:)!Gaussian Random Field
      read(27,*) tmp(6,i,:,:)! cross power spectrum of linear and non-linear delta field
      read(28,*) tmp(7,i,:,:)! cross power spectrum of linear and reconstructed delta field
      read(29,*) tmp(8,i,:,:)! cross power spectrum of linear and linear delta field @ z=0
      read(30,*) tmp(9,i,:,:)! cross power spectrum of linear and linear delta field @ z=100
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)
      close(30)
   enddo
i=0
   do xbin=1,nbin
      tmp_kbin(xbin)=tmp(1,1,1,xbin)
      if (xbin.eq.2 .or. xbin.eq.3 .or. xbin.eq.4 .or. xbin.eq.6 .or. xbin.eq.8 .or. xbin.eq.12 ) cycle
      i=i+1
      kbin(i)=tmp_kbin(xbin)
      pow(1,:,i)=tmp(1,:,2,xbin)!n-body
      pow(2,:,i)=tmp(2,:,2,xbin)!reconstruction
      pow(2,:,i)=pow(2,:,i)*kbin(i)**4
      pow(3,:,i)=tmp(3,:,2,xbin)!linear at z=0
      pow(4,:,i)=tmp(4,:,2,xbin)!linear at z=100
      pow(5,:,i)=tmp(5,:,2,xbin)!Gaussian Random Field
      pow(6,:,i)=tmp(6,:,2,xbin)!cross power spectrum of linear and non-linear delta field
      pow(7,:,i)=tmp(7,:,2,xbin)!cross power spectrum of linear and reconstructed delta field
      pow(8,:,i)=tmp(8,:,2,xbin)!cross power spectrum of linear and linear delta field @ z=0
      pow(9,:,i)=tmp(9,:,2,xbin)!cross power spectrum of linear and linear delta field @ z=100
   enddo
   do k=1,4
      pow_fine(k,:,:)=pow(k,:,:)/pow(5,:,:)
      correlation(k,:,:)=pow(k+5,:,:)/sqrt(pow(k,:,:)*pow(4,:,:))
   enddo
     correlation(5,:,:)=pow(5,:,:)/sqrt(pow(5,:,:)*pow(5,:,:))
!     correlation(1,:,:)=1
     pow_fine(5,:,:)=pow(5,:,:)/pow(5,:,:)
     
!  k=1
!  tmp_ps=1
!  do i=1,nk
!     do j=tmp_ps,nkraw
!        if ( kbin(i)>=tmp_lin(1,1,j) .and. kbin(i)<tmp_lin(1,1,j+1)) then
!           pow(3,:,k)=tmp_lin(:,2,j)*10**((log10(kbin(i))-log10(tmp_lin(1,1,j)))*(log10(tmp_lin(:,2,j+1))-log10(tmp_lin(:,2,j)))&
!                           /(log10(tmp_lin(1,1,j+1))-log10(tmp_lin(1,1,j))))
!           pow(4,:,k)=tmp_lin100(:,2,j)*10**((log10(kbin(i))-log10(tmp_lin(1,1,j)))*(log10(tmp_lin100(:,2,j+1))-log10(tmp_lin100(:,2,j)))&
!                           /(log10(tmp_lin(1,1,j+1))-log10(tmp_lin(1,1,j))))
!           k=k+1
!           tmp_ps=j
!           exit
!        else
!           cycle
!        endif
!     enddo
!     if (k==nk+1) exit
!  enddo
   
   print*,'kbins=',tmp(1,1,1,:)

   print*,'read powerspectrum done.'
!   print*,'tmp_kbins=',tmp_kbin
!   print*,'kbins=',kbin
!   print*,'power spectrum=',pow(:,1,:)
   meanpow=sum(pow(1:5,:,:),2)/n_sim
   print*,'meanpow=',meanpow
   meanpow_fine=sum(pow_fine,2)/n_sim
   meancorrelation=sum(correlation,2)/n_sim
   print*,'meanpow=',meanpow_fine
   print*,'meancorrelation=',meancorrelation
   do k=1,5
      do ybin=1,nk
      do xbin=1,nk
            pow1=pow(3,:,xbin)
            pow2=pow(3,:,ybin)
            call calccov(pow1,pow2,n_sim,cov(k,xbin,ybin),covnorm(k,xbin,ybin))
            if (xbin==ybin ) then
               if (k==5) print*,'cov=',cov(:,xbin,xbin)
               call inverse(reshape(covnorm(k,1:xbin,1:ybin),(/xbin,ybin/)),invcov(1:xbin,1:ybin),xbin)
               if (k==5) print*,'inv=',invcov(xbin,xbin)
               tmp_fishinfo=0
               do j=1,xbin
               do i=1,xbin
                      tmp_fishinfo(j)=tmp_fishinfo(j)+invcov(i,j)*meancorrelation(k,i)
               enddo
               enddo
               do j=1,xbin
                      fishinfo(k,xbin)=fishinfo(k,xbin)+tmp_fishinfo(j)*meancorrelation(k,j)
               enddo
!          
!   !            invcov=MATMUL(reshape(meanpow(k,:),(/nk,1/)),invcov)
!   !            fishinfo(k,:)=MATMUL(MATMUL(meanpow(k,:),invcov),reshape(meanpow(k,:),(/1,nk/)))
!   !            do i=1,xbin
!   !            do j=1,ybin
!   !               fishinfo(k,xbin)=fishinfo(k,xbin)+invcov(i,j)
!   !               stn(k,xbin)=sqrt(abs(stn(k,xbin)+pow1(1)*pow2(1)*invcov(i,j)))
!   !            enddo
!   !            enddo
!              ! checkinv=MATMUL(reshape(cov(1,:,:),(/nk,nk/)),reshape(invcov(1,:,:),(/xbin,ybin/)))
!   !            write(whichbin,'(i3)') xbin
!   !            print*,k,'fisherinfo(k=',trim((adjustl(whichbin))),')=',fishinfo(k,xbin)
!   !          endif
!   !     do j=1,nk-1
!   !     do i=1,nk-1
!   !            invcov(i,j)=(log(meanpow(k,i+1))-log(meanpow(k,i)))/(log(meanpow(3,i+1))-log(meanpow(3,i)))*invcov(i,j)
!   !     enddo
!   !     enddo
!   !     do j=1,nk-1
!   !     do i=1,nk-1
!   !            invcov(i,j)=invcov(i,j)*(log(meanpow(k,j+1))-log(meanpow(k,j)))/(log(meanpow(3,j+1))-log(meanpow(3,j)))
!   !     enddo
!   !     enddo
!   !     do i=1,nk
!   !            fishinfo(k,i)=invcov(i,i)
!               do i=1,xbin
!               do j=1,ybin
!                  
!                  fishinfo(k,xbin)=fishinfo(k,xbin)+invcov(i,j)
!   !               stn(k,xbin)=sqrt(abs(stn(k,xbin)+pow1(1)*pow2(1)*invcov(i,j)))
!               enddo
!               enddo
                           !               invcov(1,1:xbin)=MATMUL(reshape(meancorrelation(k,1:xbin),(/1,xbin/)),invcov(1:xbin,1:xbin))
                           !               if (k==4 .and. xbin==nk) print*, 'invcov=',invcov(1,1:xbin)
                           !               fishinfo(k,xbin)=MATMUL(invcov(1,1:xbin),reshape(meancorrelation(k,1:xbin),(/xbin,1/)))
                           !               if (k==4 .and. xbin==nk) print*, 'invcov=',invcov(1,1)
!           do i=1,xbin
!           do j=1,ybin
!               fishinfo(k,xbin)=fishinfo(k,xbin)+invcov(i,j)
!           enddo
!           enddo
               write(whichbin,'(i3)') xbin
               print*,k,'fisherinfo(k=',trim((adjustl(whichbin))),')=',fishinfo(k,xbin)
           endif
      enddo
      enddo
   enddo
   do xbin=1,nk
   write(whichbin,'(i3)') xbin
   print*,'nb,rec,lin,lin100 fisherinfo(k=',trim((adjustl(whichbin))),')=',fishinfo(:,xbin)
   !print*,'nb,rec,lin,lin100 stn(k=',trim((adjustl(whichbin))),')=',stn(:,xbin)
   enddo
   checkinv=MATMUL(reshape(covnorm(5,:,:),(/nk,nk/)),invcov)
   do ybin=1,nk
   do xbin=1,nk
      do k=1,5
         corr(k,xbin,ybin)=abs(cov(k,xbin,ybin))/sqrt(abs(cov(k,xbin,xbin)*cov(k,ybin,ybin)))
      enddo
   enddo
   enddo
   open(10,file='kbin.dat',status='replace',access='stream')
   open(11,file='meanpow_nb.dat',status='replace',access='stream')
   open(12,file='meanpow_nb_fine.dat',status='replace',access='stream')
   open(13,file='errorbar_nb_2.dat',status='replace',access='stream')
   open(14,file='meanpow_rec.dat',status='replace',access='stream')
   open(15,file='meanpow_rec_fine.dat',status='replace',access='stream')
   open(16,file='errorbar_rec_2.dat',status='replace',access='stream')
   open(17,file='meanpow_z0.dat',status='replace',access='stream')
   open(18,file='meanpow_z0_fine.dat',status='replace',access='stream')
   open(19,file='errorbar_z0_2.dat',status='replace',access='stream')
   open(20,file='meanpow_z100.dat',status='replace',access='stream')
   open(21,file='meanpow_z100_fine.dat',status='replace',access='stream')
   open(22,file='errorbar_z100_2.dat',status='replace',access='stream')
   open(23,file='meanpow_rdn.dat',status='replace',access='stream')
   open(24,file='meanpow_rdn_fine.dat',status='replace',access='stream')
   open(25,file='errorbar_rdn_2.dat',status='replace',access='stream')
do k=1,5
print*,k,'errorbar'
   do xbin=1,nk
if (k ==1 )      write(10) kbin(xbin)
         call calerror(pow(k,:,xbin),error(k,xbin),n_sim)
!         pow(k,1,xbin)=sum(pow(k,:,xbin))/n_sim
!if (k == 2)  pow(k,1,xbin)=pow(k,1,xbin)*kbin(xbin)**4
         write(3*k+8) meanpow(k,xbin)
         write(3*k+9) meanpow_fine(k,xbin)
         write(3*k+10) error(k,xbin)
         
!print*,xbin,error(k,xbin)
   enddo
enddo
   do k=1,13
      close(k+9)
   enddo
   print*,'power spectrum of nc=', pow(1,1,:)
   print*,'power spectrum of rec=', pow(2,1,:)
   print*,'power spectrum of linear=', pow(3,1,:)
   print*,'power spectrum of linear100=', pow(4,1,:)
   open(23,file='fishinfo_nb_2.dat',status='replace',access='stream')
   open(24,file='fishinfo_rec_2.dat',status='replace',access='stream')
   open(25,file='fishinfo_z0_2.dat',status='replace',access='stream')
   open(26,file='fishinfo_z100_2.dat',status='replace',access='stream')
   open(27,file='correlation_nb.dat',status='replace',access='stream')
   open(28,file='correlation_rec.dat',status='replace',access='stream')
   open(29,file='correlation_z0.dat',status='replace',access='stream')
   open(30,file='correlation_z100.dat',status='replace',access='stream')
   open(31,file='stn_nb_2.dat',status='replace',access='stream')
   open(32,file='stn_rec_2.dat',status='replace',access='stream')
   open(33,file='stn_z0_2.dat',status='replace',access='stream')
   open(34,file='stn_z100_2.dat',status='replace',access='stream')
   do k=1,4
      write(k+22) fishinfo(k,:)
      close(k+22)
   enddo
   do k=1,4
      write(k+26) corr(k,:,:)
      close(k+26)
   enddo
   do k=1,4
      write(k+30) stn(k,:)
      close(k+30)
   enddo
   print*, 'fishinfo written!'
   do i=1,nk
      print*,invcov(i,i),checkinv(i,i)
   enddo
  open(35,file='fishinfo_rdn_2.dat',status='replace',access='stream')
  open(36,file='correlation_rdn_2.dat',status='replace',access='stream')
  open(37,file='stn_rdn_2.dat',status='replace',access='stream')
  write(35) fishinfo(5,:)
  write(36) corr(5,:,:)
  write(37) stn(5,:)
  close(35)
  close(36)
  close(37)
!print*,pow(3,:,8)
end

subroutine calerror(pow,error,n_sim)
   implicit none
   integer n_sim,i
   double precision pow(n_sim),pow_tmp(n_sim),error
   error=sum(pow)/n_sim
do i=1,n_sim
   pow_tmp(i)=(pow(i)-error)**2
enddo
   error=sqrt(sum(pow_tmp)/(n_sim-1))/sqrt(real(n_sim))
endsubroutine calerror

subroutine calccov(pow1,pow2,n_sim,cov12,covnorm12)

   implicit none
   integer n_sim,i
   double precision pow1(n_sim),pow2(n_sim),cov12,covnorm12,pow11(n_sim),pow22(n_sim)
   pow11=pow1-sum(pow1)/n_sim
   pow22=pow2-sum(pow2)/n_sim
!   do i=1,n_sim
!      cov12=cov12+pow11(i)*pow22(i)
!   enddo
   cov12=DOT_PRODUCT(pow11,pow22)/(n_sim-1)
   covnorm12=cov12/sum(pow1)/sum(pow2)*n_sim**2
!   cov12=abs(cov12)
endsubroutine calccov

subroutine inverse(A,IA,N)
  implicit none
  integer N
  double precision A(N,N),IA(N,N)
  double precision,allocatable :: B(:,:)
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
!  double precision  :: matrix(:,:)
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
  double precision    :: M(N,N)
  double precision    :: S(N,N)
  integer :: I,J
  double precision :: E
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
  double precision    :: M(N,N)
  double precision    :: S(N,N)
  integer :: I,J
  double precision :: E
  do I=N,2,-1
    do J=I-1,1,-1
      E=M(J,I)/M(I,I)
      M(J,1:N)=M(J,1:N)-M(I,1:N)*E
      S(J,:)=S(J,:)-S(I,:)*E
    end do
  end do
  return
end subroutine Lower

