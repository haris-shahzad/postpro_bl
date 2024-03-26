subroutine allocate_vars()
!
! Allocate variables for the computation
!
 use mod_post
 implicit none
!
!
 allocate(dcsidx(nx),dcsidx2(nx),dcsidxs(nx))
 allocate(detady(ny),detady2(ny),detadys(ny))
 allocate(dzitdz(nz),dzitdz2(nz),dzitdzs(nz))
 allocate(x(1-ng:nx+ng))
 allocate(y(1-ng:ny+ng))
 allocate(z(1-ng:nz+ng))
 allocate(xg(1-ng:nxmax+ng+1))
 allocate(yg(1-ng:nymax+ng))
 allocate(zg(1-ng:nzmax+ng))
 allocate(xavgg(interv))
 allocate(xavg(nci))
 allocate(dU(nci))
 allocate(dxg(nxmax))
 allocate(dyg(nymax))
 allocate(dzg(nzmax))
 allocate(c(3))
 allocate(cc(0:3))
 allocate(wav(nvmean,nx,ny))
 allocate(wavg(nvmean,nci,ny))
 allocate(uy(nx,ny), uyavg(nci,ny))
 allocate(ux(nx,ny), uxavg(nci,ny))
 allocate(wavelen(int(nz/2)+1))
 allocate(wavenum(int(nz/2)+1))
 allocate(rs(6,nx,ny), rsavg(nvmean,nci,ny))
 allocate(tauwg(nxmax))
 allocate(tauw(nx), tauwavg(nci))
 allocate(rhow(nx), rhowavg(nci))
 allocate(utau(nx), utauavg(nci))
 allocate(deltav(nx), deltavavg(nci))
 allocate(cf(nx), cfavg(nci))
 allocate(d99(nx), d99avg(nci))
 allocate(retau(nx), retauavg(nci))
 allocate(dstar(nx), dstaravg(nci))
 allocate(theta(nx), thetaavg(nci))
 allocate(shapef(nx), shapefavg(nci))
 allocate(ibl(nci), ebl(nci))
!
endsubroutine allocate_vars

