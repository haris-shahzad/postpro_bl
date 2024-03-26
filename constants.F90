subroutine constants
!
! Define useful quantities
!
 use mod_post
 implicit none
!
 integer :: ii,j,jj,kk,mm
!      
 write(chx,1003) ncoords(1)
 write(chy,1003) ncoords(2)
 write(chz,1003) ncoords(3)
 1003 format(I3.3)
!
 if (visc_type==2) s2tinf = 110.4_mykind/tref_dimensional
!
 rho0 = 1._mykind
 u0   = sqrt(gamma)*rm
 v0   = 0._mykind
 w0   = 0._mykind
 s0   = 0._mykind
 p0   = 1._mykind
 t0   = 1._mykind
!
 sqgmr = sqrt(gamma)*rm/re
!
! Adiabatic wall temperature and effective wall temperature
!
 taw   = t0*(1._mykind+0.5_mykind*gm1*rm*rm*rfac)
 twall = taw*trat
!
!
! Coefficients for computation of first derivatives
!

! Evaluation of metric terms
 mm = iorder/2
 select case (mm) ! coefficient for first derivatives
 case (1)
  c(1) = 0.5_mykind
 case (2)
  c(1) = 2._mykind/3._mykind
  c(2) = -1._mykind/12._mykind
 case (3)
  c(1) = 0.75_mykind
  c(2) = -0.15_mykind
  c(3) = 1._mykind/60._mykind
 endselect

 select case (mm) ! coefficient for second derivatives
 case (1)
  cc(0) = -2._mykind
  cc(1) =  1._mykind
 case (2)
  cc(0) = -2.5_mykind
  cc(1) = 4._mykind/3._mykind
  cc(2) = -1._mykind/12._mykind
 case (3)
  cc(0) = -245._mykind/90._mykind
  cc(1) = 1.5_mykind
  cc(2) = -0.15_mykind
  cc(3) = 1._mykind/90._mykind
 endselect
!
end subroutine constants    
