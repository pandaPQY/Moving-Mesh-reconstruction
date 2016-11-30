program gridfigure

  implicit none
  integer,parameter :: nf=128
  integer,parameter :: ng=1024
  integer,parameter :: ratio=ng/nf
  real def(nf,nf,nf)
  real tmp_dsp(3,0:nf,0:nf,0:nf)
  real tmp(3,0:nf,0:nf,0:nf)
  real dsp(3,ng,ng,ng)
  real den(ng,ng,ng)
  real recden(ng,ng,ng)
  real part(2,2,2)
  integer counting(ng,ng,ng)

  integer,parameter :: nplot=1 ! thickness of the colum
  integer i,j,k,ip,jp,kp,k1,k2,im,jm,k1m,nk1,nk2
  real kk1,kk2,r1,r2
  real proj(ng,ng)
  real proj2(2,nf,nf)

  real newposition(3,ng,ng,ng)
!  integer i,j,k,ip,jp,kp,n
  real ib,jb,kb
  character(200),parameter :: inputfile1='/scratch2/p/pen/qiaoyin/datafromsimulation/fromsimulation_1024/recdef/0.000def_0128_2.bin'
  character(200),parameter :: inputfile2='/scratch2/p/pen/qiaoyin/datafromsimulation/fromsimulation_1024/0.000den1024_2.bin'
  character(200),parameter :: outputfile1='newposition_1024_proj.dat'
  character(200),parameter :: outputfile2='recden_1024_proj.dat'
  open(15,file=trim(adjustl(inputfile1)),access='stream')
  read(15) def
  print*,'reading : ',inputfile1
  close(15)
  open(15,file=trim(adjustl(inputfile2)),access='stream')
  read(15) den
  den=den/sum(den)*ng**3
  print*,'reading : ',inputfile2
  close(15)
  
  print*,'average of den = ',sum(den)/ng**3

 ! find the new position according to the deformation potential
  do k=1,nf
  do j=1,nf
  do i=1,nf
    ip=mod(i,nf)+1
    jp=mod(j,nf)+1
    kp=mod(k,nf)+1 
    tmp_dsp(1,i,j,k)=real(ratio)*1./4.*(def(ip,j,k)-def(i,j,k)+def(ip,jp,k)-def(i,jp,k)+def(ip,j,kp)-def(i,j,kp)+def(ip,jp,kp)-def(i,jp,kp))
    tmp_dsp(2,i,j,k)=real(ratio)*1./4.*(def(i,jp,k)-def(i,j,k)+def(ip,jp,k)-def(ip,j,k)+def(i,jp,kp)-def(i,j,kp)+def(ip,jp,kp)-def(ip,j,kp))
    tmp_dsp(3,i,j,k)=real(ratio)*1./4.*(def(i,j,kp)-def(i,j,k)+def(ip,j,kp)-def(ip,j,k)+def(i,jp,kp)-def(i,jp,k)+def(ip,jp,kp)-def(ip,jp,k))
!    newposition(1,8*i,8*j,8*k)=i+0.5+dsp(1,8*i,8*j,8*k)
!    newposition(2,8*i,8*j,8*k)=j+0.5+dsp(2,8*i,8*j,8*k)
!    newposition(3,8*i,8*j,8*k)=k+0.5+dsp(3,8*i,8*j,8*k)
  enddo
  enddo
  enddo
  open(11,file='dsp_proj.dat',status='replace',access='stream')
  write(11) tmp_dsp(1:2,:,:,nf/2:nf/2+3)
  close(11)
!  do i=1,3
  tmp=cshift(tmp_dsp,-1,2)
  tmp=cshift(tmp,-1,3)
  tmp=cshift(tmp,-1,4)
!  enddo
  tmp_dsp(:,0,1:nf,1:nf)=tmp(:,0,1:nf,1:nf)
  tmp_dsp(:,1:nf,0:nf,0)=tmp(:,1:nf,0:nf,0)
  tmp_dsp(:,1:nf,0,1:nf)=tmp(:,1:nf,0,1:nf)
!  tmp_dsp(0,0,0)=tmp(nf,nf,nf)

  do k=1,ng
  do j=1,ng
  do i=1,ng
     ib=mod(i,ratio)/real(ratio)
     jb=mod(j,ratio)/real(ratio)
     kb=mod(k,ratio)/real(ratio)
     ip=floor(i/real(ratio))+1-floor(ib)
     jp=floor(j/real(ratio))+1-floor(jb)
     kp=floor(k/real(ratio))+1-floor(kb)
      part(1,1,1)=ib*jb*kb
      part(2,1,1)=(1-ib)*jb*kb
      part(1,2,1)=ib*(1-jb)*kb
      part(1,1,2)=ib*jb*(1-kb)
      part(2,2,1)=(1-ib)*(1-jb)*kb
      part(2,1,2)=(1-ib)*jb*(1-kb)
      part(1,2,2)=ib*(1-jb)*(1-kb)
      part(2,2,2)=(1-ib)*(1-jb)*(1-kb)
      dsp(:,i,j,k)=dsp(:,i,j,k)+tmp_dsp(:,ip-1,jp-1,kp-1)*part(1,1,1)
      dsp(:,i,j,k)=dsp(:,i,j,k)+tmp_dsp(:,ip,jp-1,kp-1)*part(2,1,1)
      dsp(:,i,j,k)=dsp(:,i,j,k)+tmp_dsp(:,ip-1,jp,kp-1)*part(1,2,1)
      dsp(:,i,j,k)=dsp(:,i,j,k)+tmp_dsp(:,ip-1,jp-1,kp)*part(1,1,2)
      dsp(:,i,j,k)=dsp(:,i,j,k)+tmp_dsp(:,ip,jp,kp-1)*part(2,2,1)
      dsp(:,i,j,k)=dsp(:,i,j,k)+tmp_dsp(:,ip,jp-1,kp)*part(2,1,2)
      dsp(:,i,j,k)=dsp(:,i,j,k)+tmp_dsp(:,ip-1,jp,kp)*part(1,2,2)
      dsp(:,i,j,k)=dsp(:,i,j,k)+tmp_dsp(:,ip,jp,kp)*part(2,2,2)
      newposition(1,i,j,k)=i+0.5+dsp(1,i,j,k)
      newposition(2,i,j,k)=j+0.5+dsp(2,i,j,k)
      newposition(3,i,j,k)=k+0.5+dsp(3,i,j,k)
  enddo
  enddo
  enddo
  !replace the densify field to one satisfying the curcilinear grids

! get corrected column density
  k1=ng/2
  k2=k1+nplot-1
  do j=1,ng
  do i=1,ng
    jm=modulo(j-2,ng)+1
    im=modulo(i-2,ng)+1
    k1m=modulo(k1-2,ng)+1
    kk1=k1-1+(dsp(3,im,jm,k1m)+dsp(3,i,jm,k1m)+dsp(3,im,j,k1m)+dsp(3,i,j,k1m))/4
    kk2=k2+(dsp(3,im,jm,k2)+dsp(3,i,jm,k2)+dsp(3,im,j,k2)+dsp(3,i,j,k2))/4
    nk1=floor(kk1)+1
    nk2=floor(kk2)+1
    r1=nk1-kk1
    r2=kk2-nk2+1
    proj(i,j)=r1*den(i,j,nk1)+sum(den(i,j,nk1+1:nk2-1))+r2*den(i,j,nk2)-sum(den(i,j,nk2:nk1))
  enddo
  enddo
  proj=proj/nplot
  proj2=newposition(1:2,ratio::ratio,ratio::ratio,ng/2)!+newposition(1:2,ratio::ratio,ratio::ratio,ng))/3.

  open(11,file=trim(adjustl(outputfile1)),status='replace',access='stream')
  print*,'writing : ',outputfile1
  write(11) proj2
  close(11)

  open(11,file=trim(adjustl(outputfile2)),status='replace',access='stream')
  write(11) proj
  print*,'writing : ',outputfile2
  close(11)


!
!  counting=0
!  recden=0
!  n=0
!  print*,'begin to produce recden : '
!  do k=1,ng
!  do j=1,ng
!  do i=1,ng
!     ip=floor(newposition(1,i,j,k))
!     jp=floor(newposition(2,i,j,k))
!     kp=floor(newposition(3,i,j,k))
!     recden(ip,jp,k)=recden(ip,jp,k)+den(ip,jp,kp)
!     counting(ip,jp,k)=counting(ip,jp,k)+1
!  enddo
!  enddo
!  enddo
!
!  do k=1,ng
!  do j=1,ng
!  do i=1,ng
!     if ( counting(i,j,k)==0 ) counting(i,j,k)=1
!     if ( counting(i,j,k)==0 ) n=n+1
!  enddo
!  enddo
!  enddo
!  recden=recden/counting
!  print*,'max counting = ',maxval(counting),' where ? ',maxloc(counting)
!  print*,' min counting = ',minval(counting),'number of vacant = ',n
!

!  open(15,file=trim(adjustl(outputfile1)),status='replace',access='stream')
!  print*,'writing : ',outputfile1
!  write(15) newposition(:,8::8,8::8,64)
!  close(15)
!  
!  open(15,file=trim(adjustl(outputfile2)),status='replace',access='stream')
!  write(15) recden(:,:,60)
!  print*,'writing : ',outputfile2
!  close(15)  
!
!  open(15,file='dsp_1024_2.dat',status='replace',access='stream')
!  write(15) dsp
!  close(15)
!
!  print*, newposition(:,20,20,20),dsp(:,20,20,20),i,j,k,ip,jp,kp


end
