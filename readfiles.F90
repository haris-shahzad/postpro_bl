subroutine readinp
!
! Reading the input file
!
 use mod_post
 implicit none
!
 integer :: i_skip, is
!
!
 open (unit=12,file='input.dat',form='formatted')
 do i_skip=1,39
  read (12,*)
 enddo
 read (12,*)
 read (12,*) rm, retauinflow, trat, visc_type, tref_dimensional
 close(12)
!

 open (unit=12,file='reynolds.dat',form='formatted')
 read (12,*) re
 close(12)

 open (unit=12,file='fortran_input.dat',form='formatted')
 read (12,*) wavg_toggle 
 read (12,*) spectra_toggle, spec_l_toggle 
 if (spec_l_toggle == 1) then
  spectra_toggle = 0
  read(12,*) nxspec
  allocate(igslice(nxspec))
  read(12,*) (igslice(is),is=1,nxspec)
 endif
 read (12,*) nonav_toggle, imin, imax 
 read (12,*) ibl_toggle
 close(12)
!
end subroutine readinp


subroutine readgrid
!
! Reading the mesh
!
 use mod_post
 implicit none
!
 integer :: i,j,k
 integer :: i1,i2,j1,j2,k1,k2
!
! x
 i1 = 1
 i2 = nxmax+1
 open(10,file='x.dat')
 do i=i1,i2
  read(10,*) xg(i)
 enddo
 close(10)
! y
 open(10,file='y.dat')
 j1 = 1
 j2 = nymax
 do j=j1,j2
  read(10,*) yg(j)
 enddo
 close(10)
! z
 if (ndim==3) then
  k1 = 1
  k2 = nzmax
  open(10,file='z.dat')
  do k=k1,k2
   read(10,*) zg(k)
  enddo
  close(10)
 endif
!       
 do i=1,ng
  xg(1-i)     = 2._mykind*xg(2-i)-xg(3-i)
  xg(nxmax+i) = 2._mykind*xg(nxmax+i-1)-xg(nxmax+i-2)
 enddo
 xg(nxmax+ng+1) = 2._mykind*xg(nxmax+ng)-xg(nxmax+ng-1)

 do j=1,ng
  yg(1-j)     = 2._mykind*yg(2-j)-yg(3-j)
  yg(nymax+j) = 2._mykind*yg(nymax+j-1)-yg(nymax+j-2)
 enddo

 do k=1,ng
  zg(1-k)     = 2._mykind*zg(2-k)-zg(3-k)
  zg(nzmax+k) = 2._mykind*zg(nzmax+k-1)-zg(nzmax+k-2)
 enddo



!
end subroutine readgrid

subroutine readstats
 !
 use mod_post
 implicit none
!
 integer :: l,m
 integer :: mpi_io_file
 integer :: filetype
 integer :: local_datatype
 integer, dimension(3) :: sizes     ! Dimensions of the total grid
 integer, dimension(3) :: subsizes  ! Dimensions of grid local to a procs
 integer, dimension(3) :: starts    ! Starting coordinates
 integer, dimension(3) :: memsizes
 integer :: size_real
 integer :: j, ntotxy
 integer (kind=mpi_offset_kind) :: offset
 character(len=256) :: oldname, newname
!
 if (ncoords(3)==0) then
!
  sizes(1) = nblocks(1)*nx
  sizes(2) = nblocks(2)*ny
  sizes(3) = 1
  subsizes(1) = nx
  subsizes(2) = ny
  subsizes(3) = 1 
  starts(1) = 0 + ncoords(1)*subsizes(1)
  starts(2) = 0 + ncoords(2)*subsizes(2)
  starts(3) = 0
  ntotxy = nx*ny
!
  call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,iermpi)
  call mpi_type_commit(filetype,iermpi)
  call mpi_file_open(mp_cartx,'stat.bin',mpi_mode_rdonly,mpi_info_null,mpi_io_file,iermpi)
  offset = 0
  do l=1,nvmean
   call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
   call mpi_file_read_all(mpi_io_file,wav(l,1:nx,1:ny),ntotxy,mpi_prec,istatus,iermpi)
   call mpi_type_size(mpi_prec,size_real,iermpi)
   do m=1,nblocks(1)*nblocks(2)
    offset = offset+size_real*ntotxy
   enddo
  enddo
  call mpi_file_close(mpi_io_file,iermpi)
  call mpi_type_free(filetype,iermpi)


!
 endif
!
 call mpi_bcast(wav,nvmean*nx*ny,mpi_prec,0,mp_cartz,iermpi)
!
end subroutine readstats

subroutine readforce
 !
 use mod_post
 implicit none
 real(mykind), dimension(nxmax) :: forcet
 real(mykind) :: xx
 real(mykind), dimension(interv) :: forceav
 integer :: m, i, ii 
!

 if (masterproc) then
  open(12,file='forcex.dat')
  do i = 1,nxmax
   read(12,*) xx, forcet(i)
   !forcet(i) = forcet(i)/336._mykind ! Remember to change
  enddo
  close(12)

  forceav = 0._mykind
  do m = 0,interv-1
   do i = 1,nxint
    forceav(m+1) = forceav(m+1) + forcet(m*nxint+i)*dxg(m*nxint+i)
   enddo
  enddo
  forceav(1) = forceav(2)
  forceav(interv) = forceav(interv-1)



  xavgg = 0._mykind
  do m = 0,interv-1
   do i = 1,nxint
    xx = xg(m*nxint+nxint) - xg(m*nxint)
    tauwg(m*nxint+i) = forceav(m+1)/rlz/xx
    xavgg(m+1) = xavgg(m+1) + xg(m*nxint+i) 
   enddo
   xavgg(m+1) = xavgg(m+1)/nxint
  enddo

  open(12,file='force_av.dat')
  do i = 1,interv
   write(12,*) xavgg(i), forceav(i)
  enddo
  close(12)

 endif


 call mpi_bcast(tauwg,nxmax,mpi_prec,0,mpi_comm_world,iermpi)
 call mpi_bcast(xavgg,interv,mpi_prec,0,mpi_comm_world,iermpi)

 ii = nx*ncoords(1)
 do i=1,nx
  tauw(i) = tauwg(ii+i)
  rhow(i) = 1._mykind ! Remember to correct
 enddo

 ii = nci*ncoords(1)
 do i=1,nci
  xavg(i) = xavgg(ii+i)
 enddo

!
end subroutine readforce



subroutine readspectra
 !
 use mod_post
 implicit none
 character(3) :: nbx
 integer :: i, j, k, ii

 allocate(specz(nci,nvspec,ny,int(nz/2)+1))
 
 ii = nci*ncoords(1)
 do i = 1,nci
  write(nbx,1003) ii+i-1
  open(10,file='./Spectra_Averaged/speczav_'//nbx//'.bin', form='unformatted')
  read(10) specz(i,:,:,:)
  close(10)
 enddo

 do i = 1,nci
  do j = 1,ny
   do k = 2,int(nz/2)+1
    wavenum(k) = (k-1)*2*pi/rlz
    wavelen(k) = rlz/(k-1)/deltavavg(i) 
    specz(i,1,j,k) = specz(i,1,j,k)/utauavg(i)**2
    specz(i,2,j,k) = specz(i,2,j,k)/utauavg(i)**2
    specz(i,3,j,k) = specz(i,3,j,k)/utauavg(i)**2
    specz(i,4,j,k) = specz(i,4,j,k)/utauavg(i)**2
    specz(i,5,j,k) = specz(i,5,j,k)/1._mykind**2
    specz(i,6,j,k) = specz(i,6,j,k)/tauwavg(i)**2
   enddo
  enddo
 enddo

 1003  format(I3.3)
!
end subroutine readspectra


subroutine readspectra_l
 !
 use mod_post
 implicit none
 character(3) :: nbx
 integer :: i, j, k, ii
 real(mykind) :: uu, vv, ww, rho


 if (ncoords(1) == 0) then
  allocate(specz(nxspec,ny,int(nz/2)+1,nvspec))
  
  open(10,file='./gsliceyz_w.bin', access='stream', form='unformatted')
  read(10) specz(:,:,:,:)
  close(10)

  do i = 1,nxspec
   ii = igslice(i)
   do j = 1,ny
    do k = 2,int(nz/2)+1
     wavenum(k) = (k-1)*2*pi/rlz*deltavavg(ii)
     wavelen(k) = rlz/(k-1)/deltavavg(ii) 
     specz(i,j,k,1) = nzmax*specz(i,j,k,1)/1._mykind**2
     specz(i,j,k,2) = nzmax*specz(i,j,k,2)/utauavg(ii)**2
     specz(i,j,k,3) = nzmax*specz(i,j,k,3)/utauavg(ii)**2
     specz(i,j,k,4) = nzmax*specz(i,j,k,4)/utauavg(ii)**2
     specz(i,j,k,5) = nzmax*specz(i,j,k,5)/tauwavg(ii)**2
     specz(i,j,k,6) = nzmax*specz(i,j,k,6)/t0**2
    enddo
   enddo
  enddo

 endif

 1003  format(I3.3)
  100  format(20ES20.10)
!
end subroutine readspectra_l
