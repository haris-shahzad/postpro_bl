subroutine postprocess
 use mod_post
!
! Postprocessing the data
!
 if (masterproc) write(*,*) 'Reading statistics'
 call readstats
 call readforce
!
 call avgwav
 call calcbl(wavg, utauavg, tauwavg, rhowavg, deltavavg, cfavg, &
             d99avg, retauavg, dstaravg, thetaavg, shapefavg, nci)
 call calcrs

 !
 ! Uncomment lines below to calculate non-averaged utau etc
 !
  call calcbl(wav, utau, tauw, rhow, deltav, cf, &
              d99, retau, dstar, theta, shapef, nx)
 !

 call deltau
 call dudy
 call dudx
 !
 call writebl

 if (wavg_toggle==1) then
  if (masterproc) write(*,*) 'Writing statistics'
  call writestats
 endif

 if (nonav_toggle==1) then
  if (masterproc) write(*,*) 'Writing local statistics'
  call writeavfield
 endif

 if (ibl_toggle==1) then
  if (masterproc) write(*,*) 'Calculating IBL'
  call calcibl
  if (masterproc) write(*,*) 'Writing IBL'
  call writeibl
 endif

 if (spectra_toggle==1) then
  if (masterproc) write(*,*) 'Reading Spectra'
  call readspectra
  if (masterproc) write(*,*) 'Writing Spectra'
  call writespectra
 endif

 if (spec_l_toggle==1) then
  if (masterproc) write(*,*) 'Reading Spectra'
  call readspectra_l
  if (masterproc) write(*,*) 'Writing Spectra'
  call writespectra_l
 endif

!
end subroutine postprocess
