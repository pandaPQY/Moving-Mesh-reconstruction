subroutine scalar2gradient(indata,outdata,nn)
implicit none
integer(4)::nn(3)
real(4)::indata(nn(1),nn(2),nn(3))
real(4)::indata2(nn(1)+2,nn(2),nn(3))
real(4)::outdata(nn(1),nn(2),nn(3),3)
integer(4)::i,j,k
real(4)::k1,k2,k3
complex(4)::cmp
real(4),parameter::pi=atan(1.)*4.

! for kx
indata2(1:nn(1),:,:)=indata(1:nn(1),:,:)
call rlft3(indata2(1:nn(1),:,:),indata2(nn(1)+1:nn(1)+2,:,:),nn,1)
! * (A+Bi)*(ik)=(-kB+kAi)
do k=1,nn(3)
  do j=1,nn(2)
    do i=1,nn(1)+2,2
      k1=(i-1)/2
      cmp=cmplx(indata2(i,j,k),indata2(i+1,j,k))
      cmp=cmp*cmplx(0.,k1)/nn(1)*2.*pi
      indata2(i,j,k)=real(cmp)
      indata2(i+1,j,k)=aimag(cmp)
    enddo
  enddo
enddo
call rlft3(indata2(1:nn(1),:,:),indata2(nn(1)+1:nn(1)+2,:,:),nn,-1)
outdata(1:nn(1),:,:,1)=indata2(1:nn(1),:,:)*2./nn(1)/nn(2)/nn(3)

! for ky
indata2(1:nn(1),:,:)=indata(1:nn(1),:,:)
call rlft3(indata2(1:nn(1),:,:),indata2(nn(1)+1:nn(1)+2,:,:),nn,1)
! * (A+Bi)*(ik)=(-kB+kAi)
do k=1,nn(3)
  do j=1,nn(2)
    if (j.le.nn(2)/2+1) then 
      k2=j-1
    else
      k2=j-1-nn(2)
    endif
    do i=1,nn(1)+2,2
      cmp=cmplx(indata2(i,j,k),indata2(i+1,j,k))
      cmp=cmp*cmplx(0.,k2)/nn(2)*2.*pi
      indata2(i,j,k)=real(cmp)
      indata2(i+1,j,k)=aimag(cmp)
    enddo
  enddo
enddo
call rlft3(indata2(1:nn(1),:,:),indata2(nn(1)+1:nn(1)+2,:,:),nn,-1)
outdata(1:nn(1),:,:,2)=indata2(1:nn(1),:,:)*2./nn(1)/nn(2)/nn(3)

! for kz
indata2(1:nn(1),:,:)=indata(1:nn(1),:,:)
call rlft3(indata2(1:nn(1),:,:),indata2(nn(1)+1:nn(1)+2,:,:),nn,1)
! * (A+Bi)*(ik)=(-kB+kAi)
do k=1,nn(3)
  if (k.le.nn(3)/2+1) then 
    k3=k-1
  else
    k3=k-1-nn(3)
  endif
  do j=1,nn(2)
    do i=1,nn(1)+2,2
      cmp=cmplx(indata2(i,j,k),indata2(i+1,j,k))
      cmp=cmp*cmplx(0.,k3)/nn(3)*2.*pi
      indata2(i,j,k)=real(cmp)
      indata2(i+1,j,k)=aimag(cmp)
    enddo
  enddo
enddo
call rlft3(indata2(1:nn(1),:,:),indata2(nn(1)+1:nn(1)+2,:,:),nn,-1)
outdata(1:nn(1),:,:,3)=indata2(1:nn(1),:,:)*2./nn(1)/nn(2)/nn(3)

endsubroutine scalar2gradient
