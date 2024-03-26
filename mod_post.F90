module mod_post
 use mpi
 use, intrinsic :: iso_fortran_env, only : error_unit
 implicit none
 save
!
 integer, parameter :: doubtype = selected_real_kind(15,307)  ! double precision
 integer, parameter :: mykind    = doubtype
 integer, parameter :: mpi_prec = mpi_real8
 real(mykind), parameter :: tol_iter = 0.000000001_mykind
!
 integer, parameter :: nv      =  5   ! Physical variables
 integer, parameter :: nvmean  = 89
 integer, parameter :: nvspec  = 6
 integer, parameter :: nsolmax = 999999
 integer, parameter :: itmax   = 100000
 real(mykind), parameter :: pi = 4._mykind*atan(1._mykind)
!
! Averaging Parameters
!
 integer :: nxint
 integer :: interv, nci
!
! MPI related parameters 
!
 integer, parameter :: ndim = 3
 integer, parameter :: ndims = 3
 integer, dimension(:), allocatable  :: nblocks
 logical, dimension(:), allocatable  :: pbc
 integer, dimension(mpi_status_size) :: istatus
!
 integer, dimension(:), allocatable :: ncoords
 integer :: mp_cart,mp_cartx,mp_carty,mp_cartz
 integer :: mp_cartyz,mp_cartxy
 integer :: nrank,nproc,nrank_x, nrank_y, nrank_z
 integer :: ileftx,irightx,ilefty,irighty,ileftz,irightz
 integer :: iermpi, iercuda
!
 integer :: nxmax
 integer :: nymax,nymaxwr
 integer :: nzmax
 integer :: nx
 integer :: ny
 integer :: nz
 integer :: ng   ! Number of ghost nodes
 integer :: ngdf ! Number of ghost nodes for digital filtering
 integer :: io_type
 integer :: iorder
!
! Useful code variables
!
 real(mykind), parameter :: gamma  = 1.4_mykind
 real(mykind), parameter :: pr     = 0.72_mykind
 real(mykind), parameter :: gm1    = gamma-1._mykind
 real(mykind), parameter :: gm     = 1._mykind/gm1
 real(mykind), parameter :: ggmopr = gamma*gm/pr
 real(mykind), parameter :: vtexp  = 3._mykind/4._mykind
 real(mykind), parameter :: rfac   = 0.89_mykind !pr**(1._mykind/3._mykind)
 real(mykind) :: rm,re,sqgmr,retauinflow,s2tinf
 real(mykind) :: taw,trat,twall,tref_dimensional
 real(mykind) :: rtrms
 integer :: visc_type, jwall
 real(mykind) :: rho0,t0,p0,u0,v0,w0,s0
 logical :: masterproc
 logical :: dfupdated
!
 character(3) :: chx,chy,chz
 character(6) :: stat_io
!
! Coordinates and metric related quantities 
 real(mykind) :: rlx,rly,rlz
 real(mykind), dimension(:), allocatable :: x,y,z
 real(mykind), dimension(:), allocatable :: xg,yg,zg
 real(mykind), dimension(:), allocatable :: dcsidx,dcsidx2,dcsidxs
 real(mykind), dimension(:), allocatable :: detady,detady2,detadys
 real(mykind), dimension(:), allocatable :: dzitdz,dzitdz2,dzitdzs
 real(mykind), dimension(:), allocatable :: dxg,dyg,dzg
!
! Statistical quantities
 integer :: istat,itav,nstat,nstatloc
 real(mykind), dimension(:,:,:,:), allocatable :: specz
 real(mykind), dimension(:,:,:), allocatable :: wav, wavg
 real(mykind), dimension(:,:), allocatable :: uy, uyavg
 real(mykind), dimension(:,:), allocatable :: ux, uxavg
 real(mykind), dimension(:,:,:), allocatable :: rs, rsavg
 real(mykind), dimension(:), allocatable :: wavelen, wavenum
 real(mykind), dimension(:), allocatable :: xavgg, xavg
 real(mykind), dimension(:), allocatable :: dU
 real(mykind), dimension(:), allocatable :: tauwg, tauw, tauwavg
 real(mykind), dimension(:), allocatable :: rhow, rhowavg
 real(mykind), dimension(:), allocatable :: utau, utauavg
 real(mykind), dimension(:), allocatable :: deltav, deltavavg
 real(mykind), dimension(:), allocatable :: cf, cfavg
 real(mykind), dimension(:), allocatable :: d99, d99avg
 real(mykind), dimension(:), allocatable :: retau, retauavg
 real(mykind), dimension(:), allocatable :: dstar, dstaravg
 real(mykind), dimension(:), allocatable :: theta, thetaavg
 real(mykind), dimension(:), allocatable :: shapef, shapefavg
 real(mykind), dimension(:), allocatable :: ibl, ebl
!
 real(mykind), dimension(:), allocatable :: c
 real(mykind), dimension(:), allocatable :: cc
!
! Toggles
 integer :: wavg_toggle
 integer :: spectra_toggle, spec_l_toggle
 integer :: nonav_toggle, imin, imax
 integer :: ibl_toggle
!
! Spectra_l related quantities
 integer :: nxspec
 integer, dimension(:), allocatable :: igslice

end module mod_post
