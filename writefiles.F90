subroutine writestats
!
 use mod_post
 implicit none
 character(3) :: nbx
 integer :: i,j
 real(mykind) :: yd,yv,uu,vv,ww,grad,ui,un,rhon
 real(mykind) :: upup,vpvp,wpwp,upvp,upwp,vpwp
!
 if (ncoords(3) == 0) then
  do i = 1,nci
   write(nbx,1003) ncoords(1)*nci+i
   open(10,file='./Averages_2D/wavg_'//nbx//'.dat', form='formatted')
   write(10,100) xavg(i)
   !
   un   = utauavg(i)
   rhon = rhowavg(i)
   do j = 1,ny
    uu   = wavg(2,i,j)/un
    ui   = wavg(2,i,j)/u0
    vv   = wavg(3,i,j)/un
    ww   = wavg(4,i,j)/un
    grad = uyavg(i,j)/un*deltavavg(i)
    upup = rsavg(1,i,j)/un/un/rhon
    vpvp = rsavg(2,i,j)/un/un/rhon
    wpwp = rsavg(3,i,j)/un/un/rhon
    upvp = rsavg(4,i,j)/un/un/rhon
    upwp = rsavg(5,i,j)/un/un/rhon
    vpwp = rsavg(6,i,j)/un/un/rhon
    yd   = y(j)/d99avg(i)
    yv   = y(j)/deltavavg(i)
    write(10,100) yd,yv,uu,ui,vv,ww,grad,upup,vpvp,wpwp,upvp,upwp,vpwp
   enddo
   !
   close(10)
  enddo
 endif

 1003  format(I3.3)
  100  format(20ES20.10)

!
end subroutine writestats

subroutine writebl
!
 use mod_post
 implicit none
 character(3) :: nbx
 integer :: i
 real(mykind) :: dia
!
 dia = 0.08_mykind*0.5_mykind/0.375_mykind
 if (ncoords(3) == 0) then
  write(nbx,1003) ncoords(1)
  open(12,file='./BL_Properties/blavg_'//nbx//'.dat', form='formatted')
  open(14,file='./BL_Properties/dU_'//nbx//'.dat', form='formatted')
  do i = 1,nci
   write(12,100) xavg(i), d99avg(i),deltavavg(i),utauavg(i),tauwavg(i),retauavg(i), thetaavg(i), dstaravg(i)
   write(14,100) xavg(i), dia/d99avg(i), dia/deltavavg(i), dU(i)
  enddo
  close(12)
  close(14)
 endif

 1003  format(I3.3)
  100  format(20ES20.10)

!
end subroutine writebl


subroutine writeavfield
!
 use mod_post
 implicit none
 character(5) :: nbx
 integer :: i,j,m
 real(mykind) :: yd,yv,uu,vv,ww,ui,grad,un,rhon,grdx
 real(mykind) :: upup,vpvp,wpwp,upvp,upwp,vpwp
 real(mykind), dimension(nx) :: rhowl, utaul, deltavl, tauwl


 rhowl   = 0._mykind
 utaul   = 0._mykind
 tauwl   = 0._mykind
 deltavl = 0._mykind
 do m = 0,nci-1
  do i = 1,nxint
   rhowl(m*nxint+i) = rhowavg(m+1)
   utaul(m*nxint+i) = utauavg(m+1)
   tauwl(m*nxint+i) = tauwavg(m+1)
   deltavl(m*nxint+i) = deltavavg(m+1)
  enddo
 enddo
!
 if ((ncoords(1) == 0).and.(ncoords(3) == 0)) then
  do i = imin,imax
   write(nbx,1005) ncoords(1)*nx+i
   open(10,file='./Field_2D/wavg_'//nbx//'.dat', form='formatted')
   write(10,100) x(i)
   !
   un   = utaul(i)
   rhon = rhowl(i)
   do j = 1,ny
    uu   = wav(2,i,j)/un
    ui   = wav(2,i,j)/u0
    vv   = wav(3,i,j)/u0
    ww   = wav(4,i,j)/un
    grad = uy(i,j)/u0!*deltavl(i)
    grdx = ux(i,j)/u0*1._mykind!*deltavl(i)
    upup = rs(1,i,j)/un/un/rhon
    vpvp = rs(2,i,j)/un/un/rhon
    wpwp = rs(3,i,j)/un/un/rhon
    upvp = rs(4,i,j)/un/un/rhon
    upwp = rs(5,i,j)/un/un/rhon
    vpwp = rs(6,i,j)/un/un/rhon
    yd   = y(j)/d99(i)
    yv   = y(j)/deltavl(i)
    write(10,100) yd,yv,uu,ui,vv,ww,grad,grdx,upup,vpvp,wpwp,upvp,upwp,vpwp,y(j)
   enddo
   !
   close(10)
  enddo
 endif

 1005  format(I5.5)
  100  format(20ES20.10)

!
end subroutine writeavfield


subroutine writespectra
 !
 use mod_post
 implicit none
 character(3) :: nbx
 integer :: i, j, k, ii
 

 ii = nci*ncoords(1)
 do i = 1,nci
  write(nbx,1003) ii+i-1
  open(10,file='./Decomposed_Spectra/speczav_'//nbx//'.dat', form='formatted')
  do j = 1,ny
   do k = 2,int(nz/2)+1
    write(10,100) y(j)/deltavavg(i), &
                  wavelen(k), &
                  wavenum(k), &
                  wavenum(k)*specz(i,1,j,k), &
                  wavenum(k)*specz(i,2,j,k), &
                  wavenum(k)*specz(i,3,j,k), &
                  wavenum(k)*specz(i,4,j,k), &
                  wavenum(k)*specz(i,5,j,k), &
                  wavenum(k)*specz(i,6,j,k)

   enddo
  enddo
  close(10)
 enddo

 1003  format(I3.3)
  100  format(20ES20.10)
!
end subroutine writespectra

subroutine writespectra_l
 !
 use mod_post
 implicit none
 integer :: i, j, k, ii
 character(3) :: nbx
 

 if (ncoords(1) == 0) then
  do i = 1,nxspec
   write(nbx,1003) igslice(i)
   open(10,file='./Decomposed_Spectra/speczl_'//nbx//'.dat', form='formatted')
   do j = 1,ny
    do k = 2,int(nz/2)+1
     write(10,100) y(j)/deltavavg(igslice(i)), &
                   wavelen(k), &
                   wavenum(k), &
                   wavenum(k)*specz(i,j,k,1), &
                   wavenum(k)*specz(i,j,k,2), &
                   wavenum(k)*specz(i,j,k,3), &
                   wavenum(k)*specz(i,j,k,4), &
                   wavenum(k)*specz(i,j,k,5), &
                   wavenum(k)*specz(i,j,k,6)

    enddo
   enddo
   close(10)
  enddo
 endif

 1003  format(I3.3)
  100  format(20ES20.10)
!
end subroutine writespectra_l

subroutine writeibl
!
 use mod_post
 implicit none
 character(3) :: nbx
 integer :: i,j,m
!
 if (ncoords(3) == 0) then
  write(nbx,1003) ncoords(1)*nx
  open(10,file='./IBL/ibl_'//nbx//'.dat', form='formatted')
  !
  do i = 1,nci
   write(10,100) xavg(i), ibl(i), ibl(i)/d99avg(i)
  enddo
  !
  close(10)
 endif

 1003  format(I3.3)
  100  format(20ES20.10)

!
end subroutine writeibl
