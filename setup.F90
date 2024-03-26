subroutine setup
!
! Setup for the computation
!
 use mod_post
 implicit none
 logical :: file_exists
!
!===================================================
 if (masterproc) write(error_unit,*) 'Allocation of variables'
 call allocate_vars()
!===================================================
 if (masterproc) write(error_unit,*) 'Reading input'
 call readinp()
 call readgrid()
 call mpi_barrier(mpi_comm_world, iermpi)
!===================================================
 if (masterproc) write(*,*) 'Computing metrics'
 call constants()
 call computemetrics()
!===================================================
!===================================================
!
end subroutine setup
