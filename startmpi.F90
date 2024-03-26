subroutine startmpi
!
! Initialize MPI (and CUDA) environment
!
 use mod_post
 implicit none
!
 logical :: reord
 logical :: remain_dims(ndims)
 integer :: i_skip, dims(2), temp
 !integer(kind=cuda_count_kind) :: req_stacksize_gpu
!
 call mpi_init(iermpi)
 call mpi_comm_rank(mpi_comm_world,nrank,iermpi)
 call mpi_comm_size(mpi_comm_world,nproc,iermpi)

!
 allocate(ncoords(ndims))
 allocate(nblocks(ndims))
 allocate(pbc(ndims))
!
 open (unit=12,file='input.dat',form='formatted')
 do i_skip = 1,15
  read (12,*)
 enddo
 read (12,*)
 read (12,*) 
 read (12,*)
 read (12,*)
 read (12,*) rlx,rly,rlz
 read (12,*)
 read (12,*)
 read (12,*) nxmax,nymax,nzmax
 read (12,*)
 read (12,*)
 read (12,*) 
 read (12,*)
 read (12,*)
 read (12,*) ng, temp, iorder
 close(12)

 open (unit=12,file='mpidecomp.dat',form='formatted')
 read (12,*) nblocks(1), nblocks(3)
 read (12,*) nxint
 close(12)
!
 interv = int(nxmax/nxint)
 ngdf = iorder/2
!
 masterproc = .false.
 if (nrank==0) masterproc = .true.
!
 nblocks(2) = 1
!
 if(nblocks(1) <=0 .or. nblocks(3) <=0) then
     dims = [0,0]
     call MPI_Dims_create(nproc, 2, dims, iermpi)
     nblocks(1) = dims(1)
     nblocks(3) = dims(2)
     if (masterproc) write(error_unit,'(A,I0,A,I0)') 'Automatic MPI decomposition:', nblocks(1),' x ', nblocks(3)
 endif
!
 nx = nxmax/nblocks(1)
 ny = nymax/nblocks(2)
 nz = nzmax/nblocks(3)
 nci = nx/nxint

!
!
 pbc(1) = .false.
 pbc(2) = .false.
 pbc(3) = .true.
!
! Create 3D topology
!
 reord = .false.
 call mpi_cart_create(mpi_comm_world,ndims,nblocks,pbc,reord,mp_cart,iermpi)
 call mpi_cart_coords(mp_cart,nrank,ndims,ncoords,iermpi)
!
! Create 1D communicators
!
 remain_dims(1) = .true.
 remain_dims(2) = .false.
 remain_dims(3) = .false.
 call mpi_cart_sub(mp_cart,remain_dims,mp_cartx,iermpi)
 call mpi_comm_rank(mp_cartx,nrank_x,iermpi)
 call mpi_cart_shift(mp_cartx,0,1,ileftx,irightx,iermpi)
 remain_dims(2) = .true.
 remain_dims(1) = .false.
 remain_dims(3) = .false.
 call mpi_cart_sub(mp_cart,remain_dims,mp_carty,iermpi)
 call mpi_comm_rank(mp_carty,nrank_y,iermpi)
 call mpi_cart_shift(mp_carty,0,1,ilefty,irighty,iermpi)
 remain_dims(3) = .true.
 remain_dims(1) = .false.
 remain_dims(2) = .false.
 call mpi_cart_sub(mp_cart,remain_dims,mp_cartz,iermpi)
 call mpi_comm_rank(mp_cartz,nrank_z,iermpi)
 call mpi_cart_shift(mp_cartz,0,1,ileftz,irightz,iermpi)
!
 remain_dims(3) = .true.
 remain_dims(1) = .false.
 remain_dims(2) = .true.
 call mpi_cart_sub(mp_cart,remain_dims,mp_cartyz,iermpi)
!
 remain_dims(1) = .true.
 remain_dims(2) = .true.
 remain_dims(3) = .false.
 call mpi_cart_sub(mp_cart,remain_dims,mp_cartxy,iermpi)
!
end subroutine startmpi
