subroutine avgwav
!
 use mod_post

 implicit none
 integer :: i,m,l,j

 wavg = 0._mykind
 do m = 0,nci-1
  do i = 1,nxint
   do j = 1,ny
    do l = 1,nvmean
     wavg(l,m+1,j) = wavg(l,m+1,j) + wav(l,m*nxint+i,j)
    enddo
   enddo
   tauwavg(m+1) = tauw(m*nxint+i)
   rhowavg(m+1) = 1._mykind
  enddo
 enddo
 wavg = wavg/nxint


 
!
end subroutine avgwav

subroutine calcbl(wav, utau, tauw, rhow, deltav, cf, d99, retau, &
                  dstar, theta, shapef, nx)
!
 use mod_post, only : u0, visc_type, sqgmr, twall, vtexp, s2tinf, &
                      jwall, y, mykind, ny, nvmean 
 implicit none
 integer, intent(in) :: nx
 real(mykind), intent(in), dimension(nvmean,nx,ny) :: wav
 real(mykind), intent(in), dimension(nx) :: tauw, rhow 
 real(mykind), intent(inout), dimension(nx) :: utau, deltav, cf, d99, retau
 real(mykind), intent(inout), dimension(nx) :: dstar, theta, shapef
 integer :: i,j,j99
 real(mykind) :: dely
 real(mykind) :: dy,dyh,pdyn,pe
 real(mykind) :: rho,rhoe,rhop,rmuw,rnuw
 real(mykind) :: udel,uden,ue,unum,uu,uup
 real(mykind) :: sqgmr2,sqrtt,sdivt,sdivt1

!
!
! Mean boundary layer properties
!
 udel = 0.99_mykind*u0
 sqgmr2  = sqgmr*(1._mykind+s2tinf)
!
 do i=1,nx
  if (visc_type==1) then
   rmuw    = sqgmr*twall**vtexp
  else
   sqrtt   = sqrt(twall)
   sdivt   = s2tinf/twall
   sdivt1  = 1._mykind+sdivt
   rmuw    = sqgmr2*sqrtt/sdivt1  ! molecular viscosity
  endif
!
  utau(i)   = sqrt(abs(tauw(i))/rhow(i))
  rnuw      = rmuw/rhow(i)
  deltav(i) = rnuw/utau(i)
  pdyn      = 0.5_mykind*u0**2
  cf(i)     = tauw(i)/pdyn
  j99       = 1

  do j=1,ny-1
   uu = wav(2,i,j)
   if (uu>udel) then
    j99 = j-1
    exit
   endif
  enddo

  dely     = y(j99+1)-y(j99)
  unum     = udel-wav(2,i,j99)
  uden     = wav(2,i,j99+1)-wav(2,i,j99)
  d99(i)   = y(j99)+dely*(unum/uden) ! bl thickness
  retau(i) = d99(i)/deltav(i)
!
! Integral boundary layer thicknesses
!
  dstar(i) = 0.0_mykind
  theta(i) = 0.0_mykind
  rhoe     = wav(1,i,j99)
  pe       = wav(5,i,j99)
  ue       = wav(2,i,j99)

  do j=jwall,j99
   rho  = wav(1,i,j)/rhoe
   rhop = wav(1,i,j+1)/rhoe
   uu   = wav(2,i,j)/ue
   uup  = wav(2,i,j+1)/ue
   dy   = y(j+1)-y(j)
   dyh  = 0.5_mykind*dy
!  Trapezoidal rule
   dstar(i) = dstar(i) + dyh*((1._mykind-rho*uu)+(1._mykind-rhop*uup))
   theta(i) = theta(i) + dyh*((rho*uu*(1._mykind-uu))+(rhop*uup*(1._mykind-uup)))
  enddo
  shapef(i) = dstar(i)/theta(i) ! Shape factor H

 enddo
!
end subroutine calcbl



subroutine calcrs
!
 use mod_post
!
 implicit none
 integer :: i,j
!
 do i = 1,nx
  do j = 1,ny
   rs(1,i,j) = wav(16,i,j) - wav(1,i,j)*wav(2,i,j)*wav(2,i,j) ! rhouu
   rs(2,i,j) = wav(17,i,j) - wav(1,i,j)*wav(3,i,j)*wav(3,i,j) ! rhovv
   rs(3,i,j) = wav(18,i,j) - wav(1,i,j)*wav(4,i,j)*wav(4,i,j) ! rhoww
   rs(4,i,j) = wav(19,i,j) - wav(1,i,j)*wav(2,i,j)*wav(3,i,j) ! rhouv
   rs(5,i,j) = wav(21,i,j) - wav(1,i,j)*wav(2,i,j)*wav(4,i,j) ! rhouw
   rs(6,i,j) = wav(22,i,j) - wav(1,i,j)*wav(3,i,j)*wav(4,i,j) ! rhovw
  enddo
 enddo


 do i = 1,nci
  do j = 1,ny
   rsavg(1,i,j) = wavg(16,i,j) - wavg(1,i,j)*wavg(2,i,j)*wavg(2,i,j) ! rhouu
   rsavg(2,i,j) = wavg(17,i,j) - wavg(1,i,j)*wavg(3,i,j)*wavg(3,i,j) ! rhovv
   rsavg(3,i,j) = wavg(18,i,j) - wavg(1,i,j)*wavg(4,i,j)*wavg(4,i,j) ! rhoww
   rsavg(4,i,j) = wavg(19,i,j) - wavg(1,i,j)*wavg(2,i,j)*wavg(3,i,j) ! rhouv
   rsavg(5,i,j) = wavg(21,i,j) - wavg(1,i,j)*wavg(2,i,j)*wavg(4,i,j) ! rhouw
   rsavg(6,i,j) = wavg(22,i,j) - wavg(1,i,j)*wavg(3,i,j)*wavg(4,i,j) ! rhovw
  enddo
 enddo

!
end subroutine calcrs

subroutine deltau
!
 use mod_post
!
 implicit none
 real(mykind), dimension(1:314) :: yvelc, uvelc
 real(mykind), dimension(1:ny) :: yvel, uvel
 real(mykind) :: ys, us
 real(mykind) :: yr, ur
 integer :: i,j,ii

 if (ncoords(3) == 0) then
  yvel = 0._mykind
  uvel = 0._mykind
  ys   = 100._mykind
  yr   = 100._mykind

  open(10,file='./uvel.dat', form='formatted')
  do j = 1,314
   read(10,*) yvelc(j), uvelc(j)
  enddo
  close(10)

  call locateval(yvelc,314,ys,ii)

  us = (ys-yvelc(ii))/(yvelc(ii+1)-yvelc(ii))
  us = us*(uvelc(ii+1)-uvelc(ii)) + uvelc(ii)


  do i = 1,nci
   do j = 1,ny
    yvel(j) = y(j)/deltavavg(i)
    uvel(j) = wavg(2,i,j)/utauavg(i)
   enddo

   call locateval(yvel,ny,yr,ii)

   ur = (yr-yvel(ii))/(yvel(ii+1)-yvel(ii))
   ur = ur*(uvel(ii+1)-uvel(ii)) + uvel(ii)

   dU(i) = us - ur
  enddo
 endif
!
end subroutine deltau


subroutine dudx
!
 use mod_post
!
 implicit none
 integer :: i,j,m,l
 real(mykind) :: ccl

 if (ncoords(3) == 0) then
  ux    = 0._mykind
  uxavg = 0._mykind

  do j = 1,ny
   do i = 1,nx
    !
    do l=1,iorder/2
     ccl = c(l)
     ux(i,j) = uxavg(i,j)+ccl*(wav(2,i+l,j)-wav(2,i-l,j))
    enddo
    ux(i,j) = uxavg(i,j)*dcsidx(i)
    !
   enddo
  enddo

 endif
!
end subroutine dudx

subroutine dudy
!
 use mod_post
!
 implicit none
 integer :: i,j,m,l
 real(mykind) :: ccl
 real(mykind), dimension(1-ng:ny+ng) :: uvel

 if (ncoords(3) == 0) then
  uy    = 0._mykind
  uyavg = 0._mykind

  do i = 1,nci

   uvel(1:ny) = wavg(2,i,:)
   do m = 1,ng
    uvel(1-ng)  = -uvel(ng)
    uvel(ny+ng) =  uvel(ny)
   enddo

   do j = 1,ny
    !
    do l=1,iorder/2
     ccl = c(l)
     uyavg(i,j) = uyavg(i,j)+ccl*(uvel(j+l)-uvel(j-l))
    enddo
    uyavg(i,j) = uyavg(i,j)*detady(j)
    !
   enddo
  enddo




  do i = 1,nx

   uvel(1:ny) = wav(2,i,:)
   do m = 1,ng
    uvel(1-ng)  = -uvel(ng)
    uvel(ny+ng) =  uvel(ny)
   enddo

   do j = 1,ny
    !
    do l=1,iorder/2
     ccl = c(l)
     uy(i,j) = uy(i,j)+ccl*(uvel(j+l)-uvel(j-l))
    enddo
    uy(i,j) = uy(i,j)*detady(j)
    !
   enddo
  enddo

 endif
!
end subroutine dudy


subroutine calcibl
!
 use mod_post
!
 implicit none
 integer :: i,j,jj,j99,m,l
 integer, parameter :: s1 = 55
 integer, parameter :: s2 = 99 
 integer, parameter :: us = 48
 integer, parameter :: ds = 165
 real(mykind), dimension(ny)  :: uvelus, uvelds, uvel
 real(mykind), dimension(ny)  :: yvelus, yvelds, yvel
 real(mykind), dimension(ny)  :: uveliu, uvelid
 real(mykind) :: iblt, eblt
 real(mykind) :: uveltu, yveltu
 real(mykind) :: uveltd, yveltd
 real(mykind) :: uu, udel, dely, unum, uden

 if (ncoords(3) == 0) then
  uvelus = wavg(2,us,:)/u0
  uvelds = wavg(2,ds,:)/u0


  do i = 1,nci
   uvel = wavg(2,i,:)/u0 

   iblt = 0._mykind
   eblt = 0._mykind


   do j=jwall,ny-1
    uu = uvel(j)
    !if ((uu>0.99*uvelus(j)).and.(y(j)>0.035_mykind)) then
    if ((uu>0.9875*uvelus(j)).and.(y(j)>0.035_mykind)) then
     j99 = j-1
     exit
    endif
   enddo


   udel = 0.9875*uvelus(j99)!0.99*uvelus(j99)

   dely     = y(j99+1)-y(j99)
   unum     = udel-uvel(j99)
   uden     = uvel(j99+1)-uvel(j99)
   ibl(i)   = y(j99)+dely*(unum/uden) ! Internal Boundary Layer Thickness
   ebl(i)   = 0._mykind


  enddo

 endif

!
end subroutine calcibl




subroutine calcibl_interp
!
 use mod_post
!
 implicit none
 integer :: i,j,jj,j99,m,l
 integer, parameter :: s1 = 160
 integer, parameter :: s2 = 161
 integer, parameter :: us = 48
 integer, parameter :: ds = 165
 real(mykind), dimension(ny)  :: uvelus, uvelds, uvel
 real(mykind), dimension(ny)  :: yvelus, yvelds, yvel
 real(mykind), dimension(ny)  :: uveliu, uvelid
 real(mykind) :: iblt, eblt
 real(mykind) :: uveltu, yveltu
 real(mykind) :: uveltd, yveltd
 real(mykind) :: uu, udel, dely, unum, uden

 if (ncoords(3) == 0) then
  uvelus = wavg(2,us,:)/u0
  uvelds = wavg(2,ds,:)/u0

  yvelus = y/d99avg(us)
  yvelds = y/d99avg(ds)

  do i = s1,s2
   uvel = wavg(2,i,:)/u0 
   yvel = y/d99avg(i) 

   iblt = 0._mykind
   eblt = 0._mykind

   uveliu = 0._mykind
   uvelid = 0._mykind

   do j = jwall,ny-5
    yveltu = yvelus(j)
      
    call locateval(yvel,ny,yveltu,jj)

    uveliu(j) = (yveltu-yvel(jj))/(yvel(jj+1)-yvel(jj))
    uveliu(j) = uveliu(j)*(uvel(jj+1)-uvel(jj)) + uvel(jj)

    yveltd = yvelds(j)
      
    call locateval(yvel,ny,yveltd,jj)

    uvelid(j) = (yveltd-yvel(jj))/(yvel(jj+1)-yvel(jj))
    uvelid(j) = uvelid(j)*(uvel(jj+1)-uvel(jj)) + uvel(jj)


   enddo


   do j=jwall,ny-1
    uu = uveliu(j)
    if ((uu>0.99*uvelus(j)).and.(yvelus(j)>0.005_mykind)) then
     j99 = j-1
     print*, yvelus(j99), uveliu(j99), 0.99*uvelus(j99)
     exit
    endif
   enddo

   udel = 0.99*uvelus(j99)

   dely     = yvelus(j99+1)-yvel(j99)
   unum     = udel-uveliu(j99)
   uden     = uveliu(j99+1)-uveliu(j99)
   ibl(i)   = yvelus(j99)+dely*(unum/uden) ! Internal Boundary Layer Thickness

   do j=jwall,ny-1
    uu = uvelid(j)
    if ((uu>1.01*uvelds(j)).and.(yvelds(j)>0.005_mykind)) then
     j99 = j-1
     exit
    endif
   enddo

   udel = 1.01*uvelds(j99)

   dely     = yvelds(j99+1)-yvel(j99)
   unum     = udel-uvelid(j99)
   uden     = uvelid(j99+1)-uvelid(j99)
   ebl(i)   = yvelds(j99)+dely*(unum/uden) ! Internal Boundary Layer Thickness
   !print*, ebl(i), ibl(i)

  enddo

 endif

!
end subroutine calcibl_interp
