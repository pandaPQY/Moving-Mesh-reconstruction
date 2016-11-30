subroutine defp2delta(indata,outdata,nn)
implicit none
integer(4)::nn(3)
real(4)::indata(nn(1),nn(2),nn(3))
real(4)::indata2(nn(1)+2,nn(2),nn(3))
real(4)::outdata(nn(1),nn(2),nn(3))
integer(4)::i,j,k
real(4)::k1,k2,k3
real(4),parameter::pi=atan(1.)*4.

indata2(1:nn(1),:,:)=indata(1:nn(1),:,:)
indata2(nn(1)+1:nn(1)+2,:,:)=0.
call rlft3(indata2(1:nn(1),:,:),indata2(nn(1)+1:nn(1)+2,:,:),nn,1)
do k=1,nn(3)
  if (k.le.nn(3)/2+1) then
    k3=k-1
  else
    k3=k-1-nn(3)
  endif
  k3=k3/nn(3)*2.*pi
  do j=1,nn(2)
    if (j.le.nn(2)/2+1) then
      k2=j-1
    else
      k2=j-1-nn(2)
    endif
    k2=k2/nn(2)*2.*pi
    do i=1,nn(1)+2,2
      k1=(i-1)/2
      k1=k1/nn(1)*2.*pi

      indata2(i:i+1,j,k)=indata2(i:i+1,j,k)*(k1**2+k2**2+k3**2)

    enddo
  enddo
enddo

call rlft3(indata2(1:nn(1),:,:),indata2(nn(1)+1:nn(1)+2,:,:),nn,-1)
outdata(1:nn(1),:,:)=indata2(1:nn(1),:,:)*2./nn(1)/nn(2)/nn(3)

endsubroutine defp2delta
