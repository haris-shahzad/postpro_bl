subroutine computemetrics
!
! Computing mesh metrics
!
 use mod_post
 implicit none
!
 integer, dimension(6) :: itag
!
 real(mykind), dimension(nxmax) :: d2xg
 real(mykind), dimension(nymax) :: d2yg
 real(mykind), dimension(nzmax) :: d2zg
 real(mykind), dimension(0:nxmax) :: dxhg
 real(mykind), dimension(0:nymax) :: dyhg
 real(mykind), dimension(0:nzmax) :: dzhg
!
 real(mykind), dimension(4)   :: cs
!
 integer :: i,ii,j,jj,k,kk,l,mm,iend,jend,kend
 real(mykind) :: xloc
!
! z
 if (ndim==3) then
 else
  zg   = 0._mykind
  dzg  = 0._mykind
  d2zg = 0._mykind
 endif
!
 rlx = xg(nxmax)-xg(1)
 rly = yg(nymax)-yg(1)
 rlz = zg(nzmax+1)-zg(1)
!
! local coordinates (nodes)
 ii = nx*ncoords(1)
 do i=1-ng,nx+ng
  x(i) = xg(ii+i)
 enddo
 jj = ny*ncoords(2)
 do j=1-ng,ny+ng
  y(j) = yg(jj+j)
 enddo
 kk = nz*ncoords(3)
 do k=1-ng,nz+ng
  z(k) = zg(kk+k)
 enddo
!
 mm = iorder/2

 dxg = 0._mykind
 do i=1,nxmax
  do l=1,mm
   dxg(i) = dxg(i)+c(l)*(xg(i+l)-xg(i-l))
  enddo
  d2xg(i) = cc(0)*xg(i)
  do l=1,mm
   d2xg(i) = d2xg(i)+cc(l)*(xg(i+l)+xg(i-l))
  enddo
 enddo

 dyg = 0._mykind
 do j=1,nymax
  do l=1,mm
   dyg(j) = dyg(j)+c(l)*(yg(j+l)-yg(j-l))
  enddo
  d2yg(j) = cc(0)*yg(j)
  do l=1,mm
   d2yg(j) = d2yg(j)+cc(l)*(yg(j+l)+yg(j-l))
  enddo
 enddo

 if (ndim==3) then
  dzg = 0._mykind
  do k=1,nzmax
   do l=1,mm
    dzg(k) = dzg(k)+c(l)*(zg(k+l)-zg(k-l))
   enddo
   d2zg(k) = cc(0)*zg(k)
   do l=1,mm
    d2zg(k) = d2zg(k)+cc(l)*(zg(k+l)+zg(k-l))
   enddo
  enddo
 endif
!
 ii = nx*ncoords(1)
 do i=1,nx
  dcsidx (i) = 1._mykind/(dxg(ii+i))
  dcsidxs(i) = dcsidx(i)*dcsidx(i)
  dcsidx2(i) = -d2xg(ii+i)*dcsidxs(i)
 enddo
 jj = ny*ncoords(2)
 do j=1,ny
  detady (j) = 1._mykind/dyg(jj+j)
  detadys(j) = detady(j)*detady(j)
  detady2(j) = -d2yg(jj+j)*detadys(j)
 enddo
 if (ndim==3) then
  kk = nz*ncoords(3)
  do k=1,nz
   dzitdz (k) = 1._mykind/dzg(kk+k)
   dzitdzs(k) = dzitdz(k)*dzitdz(k)
   dzitdz2(k) = -d2zg(kk+k)*dzitdzs(k)
  enddo
 else
  do k=1,nz
   dzitdz (k) = 0._mykind
   dzitdzs(k) = 0._mykind
   dzitdz2(k) = 0._mykind
  enddo
 endif
!
!
 if (masterproc) then
  open(18,file='dxg.dat')
  do i=1,nxmax
   write(18,*) xg(i),dxg(i)
  enddo
  close(18)
  open(18,file='dyg.dat')
  do j=1,nymax
   write(18,*) yg(j),dyg(j)
  enddo
  close(18)
  if (ndim==3) then
   open(18,file='dzg.dat')
   do k=1,nzmax
    write(18,*) zg(k),dzg(k)
   enddo
   close(18)
  endif
 endif

 jwall = 1
 do j=1,ny-1
  if ((y(j)<0._mykind).and.(y(j+1)>0._mykind)) jwall = j+1
 enddo

!
end subroutine computemetrics
