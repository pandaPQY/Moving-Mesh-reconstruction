program den_assignment
implicit none
!integer,parameter :: np=223479 ! total number of particles
!integer,parameter :: nf=4608 ! number of fine cells per dim
!integer,parameter :: ng=128 ! number of grids per dim in output resolution
integer,parameter :: nn=2
integer,parameter :: nnp=256 ! particle per node per dimension
integer,parameter :: np=(nnp*nn)**3 ! total number of particles
integer,parameter :: nd=512 !distance per node
integer,parameter :: nf=nd*nn ! number of fine cells per dim
integer,parameter :: ng=1024 ! number of grids per dim in output resolution
integer n_p(nn**3) ! number of particle in one node
integer n_pt ! totally number of particle read in
real garbage(11,nn,nn,nn)
!character(len=200),parameter :: inputfile='/scratch2/p/pen/qiaoyin/simulation/simulation_2/node0/0.000xv0.dat'
character(len=200) tmp_fn
character(len=100),parameter :: dir='/scratch2/p/pen/qiaoyin/simulation/simu_'
integer,parameter :: sim=5
character(len=3) i_sim
character(len=10),parameter :: str='/node'
character(len=3) rank
character(len=20),parameter :: fn='/0.000xv'
character(len=10),parameter :: prefix='.dat'
character(len=20),parameter :: fn2='0.000den'
character(len=1),parameter :: a='_'
character(len=10),parameter :: prefix2='.bin'
character(len=4) gstr
integer i,j,k,ip,ig,pp,inn
integer i_s,j_s,k_s,ii_s,jj_s,kk_s,r,n
integer wp(np),pp_s
integer hoc_g(ng,ng,ng),ll_p(np) ! linked list
integer hoc_p(np),ll_g(ng**3) ! inverse linked list
real xv(6,np)
real tmp0(6,np)
real tmp1(6,np)
real tmp2(6,np)
real tmp3(6,np)
real tmp4(6,np)
real tmp5(6,np)
real tmp6(6,np)
real tmp7(6,np)
real xp(3,np),den(ng,ng,ng),den_ngp(ng,ng,ng)
real r2,r2min,gpos(3),hpos(3),dpos(3)

logical search
print*,'nf = ',nf
print*,'ng = ',ng
write(gstr,'(i4.4)') ng
! read in particle positions
!open(11,file='xhalo.dat',status='old',access='stream')
!read(11) xp
!close(11)
inn=0
n_pt=0
do k=1,nn
do j=1,nn
do i=1,nn
   write(rank,'(i3)') inn
   write(i_sim,'(i3)') sim
   print*,'reading node ',inn
   tmp_fn=trim(adjustl(dir))//trim(adjustl(i_sim))//trim(adjustl(str))//trim(adjustl(rank))//&
          trim(adjustl(fn))//trim(adjustl(rank))//trim(adjustl(prefix))
   print*,'reading',tmp_fn
   open(11,file=trim(adjustl(tmp_fn)),status='old',access='stream')
   read(11) n_p(inn+1)
   n_pt=n_pt+n_p(inn+1)
   read(11) garbage(:,i,j,k)
   read(11) xv(:,sum(n_p(1:inn))+1:sum(n_p(1:inn+1)))
   xv(1,sum(n_p(1:inn))+1:sum(n_p(1:inn+1)))=xv(1,sum(n_p(1:inn))+1:sum(n_p(1:inn+1)))+(i-1)*nd
   xv(2,sum(n_p(1:inn))+1:sum(n_p(1:inn+1)))=xv(2,sum(n_p(1:inn))+1:sum(n_p(1:inn+1)))+(j-1)*nd
   xv(3,sum(n_p(1:inn))+1:sum(n_p(1:inn+1)))=xv(3,sum(n_p(1:inn))+1:sum(n_p(1:inn+1)))+(k-1)*nd
   close(11)
   inn=inn+1
enddo
enddo
enddo
xp=xv(1:3,:)
print*,'max position',maxval(xp),'min posiion',minval(xp)
do ip=1,np
  xp(1,ip)=modulo(xp(1,ip),real(nf))
  xp(2,ip)=modulo(xp(2,ip),real(nf))
  xp(3,ip)=modulo(xp(3,ip),real(nf))
enddo
print*,'real position done'
print*,'max position',maxval(xp),'min posiion',minval(xp)
xp=xp*ng/nf
wp=1 ! weight -- number of grids assigned to the particle

! assign halos to grids
do ip=1,np
  i=floor(xp(1,ip))+1
  j=floor(xp(2,ip))+1
  k=floor(xp(3,ip))+1
!  if ( i<0 .or. i>ng ) print*,i
!  if ( j<0 .or. j>ng ) print*,j
!  if ( k<0 .or. k>ng ) print*,k
  den_ngp(i,j,k)=den(i,j,k)+1
  ll_p(ip)=hoc_g(i,j,k)
  hoc_g(i,j,k)=ip
enddo
print*,'rescale done'
print*,'max position',maxval(xp),'min posiion',minval(xp)

! assign non-empty grids to halos
do k=1,ng
do j=1,ng
do i=1,ng
  ig=ng**2*(k-1)+ng*(j-1)+i
  pp=hoc_g(i,j,k)
  do
    if (pp==0) exit
    hoc_p(pp)=ig
    pp=ll_p(pp)
  enddo
enddo
enddo
enddo

! assign empty grids to halos
do k=1,ng
   print*,'k = ',k
do j=1,ng
do i=1,ng
  gpos=(/i,j,k/)-0.5
  ig=ng**2*(k-1)+ng*(j-1)+i
  if (hoc_g(i,j,k)==0) then
    r=0
    search=.true.
    r2min=1000*ng**2
    do ! search larger layer
      if (.not.search) exit
      r=r+1
      do k_s=k-r,k+r
      do j_s=j-r,j+r
      do i_s=i-r,i+r
        kk_s=modulo(k_s-1,ng)+1
        jj_s=modulo(j_s-1,ng)+1
        ii_s=modulo(i_s-1,ng)+1
        pp=hoc_g(ii_s,jj_s,kk_s)
        n=0
        do
          if (pp==0) exit ! find next cell
          n=n+1
          if (n ==1000) print*,'attension',i,j,k,pp
          search=.false. ! found particle
          hpos=xp(:,pp)
          dpos=hpos-gpos
          dpos=modulo(dpos+ng/2,real(ng))-ng/2
          r2=sum(dpos**2)
          if (r2<r2min) then
            r2min=r2
            pp_s=pp
          endif
          pp=ll_p(pp)
        enddo
      enddo
      enddo
      enddo
    enddo ! finished searching current layer
    ! found nearest particle, exit searching loop
    ll_g(ig)=hoc_p(pp_s) ! let current grid point to pp_s's chain
    hoc_p(pp_s)=ig ! replace pp_s's hoc to be current grid
    wp(pp_s)=wp(pp_s)+1
  endif
enddo
enddo
enddo

! assign density
do ip=1,np
  ig=hoc_p(ip)
  do
    if (ig==0) exit
    k=(ig-1)/ng**2+1
    j=(ig-1-ng**2*(k-1))/ng+1
    i=modulo(ig-1,ng)+1
    den(i,j,k)=den(i,j,k)+1./real(wp(ip));
    ig=ll_g(ig);
  enddo
enddo

print*, 'np =',np
print*, 'sum(den) =',sum(den*1.d0)
print*, 'minval(den) =',minval(den)
den=real(den,8)/sum(real(den,8))*float(ng)**3
print*, 'np =',np
print*, 'sum(den) =',sum(den*1.d0)
print*, 'average(den) =',sum(den*1.d0)/float(ng)**3
print*, 'minval(den) =',minval(den)
tmp_fn=trim(adjustl(fn2))//trim(adjustl(gstr))//trim(adjustl(a))//trim(adjustl(i_sim))//trim(adjustl(prefix2))
open(11,file=trim(adjustl(tmp_fn)),status='replace',access='stream')
print*,'writing',tmp_fn
write(11) den
close(11)
print*, 'wrote into file den_interp.dat'

end
