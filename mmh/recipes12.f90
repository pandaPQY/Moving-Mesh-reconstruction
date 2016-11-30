SUBROUTINE fourn(data5,nn,ndim,isign)

INTEGER isign,ndim,nn(ndim)
REAL data5(*)
INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
REAL tempi,tempr
DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
ntot=1
do idim=1,ndim
  ntot=ntot*nn(idim)
end do
nprev=1
do idim=1,ndim
  n=nn(idim)
  nrem=ntot/(n*nprev)
  ip1=2*nprev
  ip2=ip1*n
  ip3=ip2*nrem
  i2rev=1
  do i2=1,ip2,ip1
    if(i2.lt.i2rev)then
      do i1=i2,i2+ip1-2,2
        do i3=i1,ip3,ip2
          i3rev=i2rev+i3-i2
          tempr=data5(i3)
          tempi=data5(i3+1)
          data5(i3)=data5(i3rev)
          data5(i3+1)=data5(i3rev+1)
          data5(i3rev)=tempr
          data5(i3rev+1)=tempi
        end do
      end do
    endif
    ibit=ip2/2
    do while( ((ibit.ge.ip1).and.(i2rev.gt.ibit)) )
      i2rev=i2rev-ibit
      ibit=ibit/2
    end do
    i2rev=i2rev+ibit
  end do
  ifp1=ip1
  do while(ifp1.lt.ip2)
    ifp2=2*ifp1
    theta=isign*6.28318530717959d0/(ifp2/ip1)
    wpr=-2.d0*sin(0.5d0*theta)**2
    wpi=sin(theta)
    wr=1.d0
    wi=0.d0
    do i3=1,ifp1,ip1
      do i1=i3,i3+ip1-2,2
        do i2=i1,ip3,ifp2
          k1=i2
          k2=k1+ifp1
          tempr=sngl(wr)*data5(k2)-sngl(wi)*data5(k2+1)
          tempi=sngl(wr)*data5(k2+1)+sngl(wi)*data5(k2)
          data5(k2)=data5(k1)-tempr
          data5(k2+1)=data5(k1+1)-tempi
          data5(k1)=data5(k1)+tempr
          data5(k1+1)=data5(k1+1)+tempi
        end do
      end do
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
    end do
    ifp1=ifp2
  end do
  nprev=n*nprev
end do

return
END

!------------------------------------------


SUBROUTINE rlft3(rdata,speq,nn,isign)

INTEGER isign,nn(3)
COMPLEX rdata(nn(1)/2,nn(2),nn(3)),speq(nn(2),nn(3))

INTEGER i1,i2,i3,j1,j2,j3
DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
COMPLEX c1,c2,h1,h2,w
c1=cmplx(0.5,0.0)
c2=cmplx(0.0,-0.5*isign)
theta=6.28318530717959d0/dble(isign*nn(1))
wpr=-2.0d0*dsin(0.5d0*theta)**2
wpi=dsin(theta)
if(isign.eq.1)then
  nn(1)=nn(1)/2
  call fourn(rdata,nn,3,isign)
  nn(1)=nn(1)*2
  do i3=1,nn(3)
    do i2=1,nn(2)
      speq(i2,i3)=rdata(1,i2,i3)
    end do
  end do
endif
do i3=1,nn(3)
  j3=1
  if (i3.ne.1) j3=nn(3)-i3+2
  wr=1.0d0
  wi=0.0d0
  do i1=1,nn(1)/4+1
    j1=nn(1)/2-i1+2
    do i2=1,nn(2)
      j2=1
      if (i2.ne.1) j2=nn(2)-i2+2
      if (i1.eq.1)then
        h1=c1*(rdata(1,i2,i3)+conjg(speq(j2,j3)))
        h2=c2*(rdata(1,i2,i3)-conjg(speq(j2,j3)))
        rdata(1,i2,i3)=h1+h2
        speq(j2,j3)=conjg(h1-h2)
      else
        h1=c1*(rdata(i1,i2,i3)+conjg(rdata(j1,j2,j3)))
        h2=c2*(rdata(i1,i2,i3)-conjg(rdata(j1,j2,j3)))
        rdata(i1,i2,i3)=h1+w*h2
        rdata(j1,j2,j3)=conjg(h1-w*h2)
      endif
    end do
    wtemp=wr
    wr=wr*wpr-wi*wpi+wr
    wi=wi*wpr+wtemp*wpi+wi
    w=cmplx(sngl(wr),sngl(wi))
  end do
end do
if(isign.eq.-1)then
  nn(1)=nn(1)/2
  call fourn(rdata,nn,3,isign)
  nn(1)=nn(1)*2
end if

return
END


