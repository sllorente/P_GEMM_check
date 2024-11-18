program test_pzgemm
   ! Given a N2xN2 distributed matrix B and a (N1+N2)x(N1+N2) distributed matrix C
   ! where:
   ! C(:,:) = 1
   ! B(:,:) = N1
   !
   ! Calculate:
   ! B(1:N2, 1:N2) = C(N1+1:N1+N2, 1:N1) * C(1:N1, N1+1:N1+N2) - B(1:N2, 1:N2)
   !
   ! The expected result is B(:,:) = 0
   ! The program check the result and print the elements of B if they are not 0.

   use mpi_f08
   use iso_fortran_env, only: sp=>real32, dp=>real64, stdout=>output_unit, int32 
   use scalapack_interface

   implicit none

   integer, parameter :: mb=4, nb=4, N1=4*5, N2=2**10 + 1 
   integer, parameter :: MAX_PRINTED_LINES=20

#if defined(fD)
#define PGEMM pdgemm
   real(kind=dp), allocatable :: C(:), B(:)
   real(kind=dp) :: alpha, beta
   real(kind=dp), parameter :: ONE = 1.0_dp, ZERO = 0.0_dp
#elif defined(fC) 
#define PGEMM pcgemm
   complex(kind=sp), allocatable :: C(:), B(:)
   complex(kind=sp) :: alpha, beta
   complex(kind=sp), parameter :: ONE = (1.0_sp, 0.0_sp), ZERO = (0.0_sp, 0.0_sp)
#elif defined(fZ)   
#define PGEMM pzgemm
   complex(kind=dp), allocatable :: C(:), B(:)
   complex(kind=dp) :: alpha, beta
   complex(kind=dp), parameter :: ONE = (1.0_dp, 0.0_dp), ZERO = (0.0_dp, 0.0_dp)
#else
#define PGEMM psgemm
   real(kind=sp), allocatable :: C(:), B(:)
   real(kind=sp) :: alpha, beta
   real(kind=sp), parameter :: ONE = 1.0_sp, ZERO = 0.0_sp
#endif

   integer :: Cdesc(dlen_), Bdesc(dlen_)
   integer :: ictxt, nprow, npcol,  myprow, mypcol

   integer(kind=int32) :: np, ierror, myrank, dims(2)

   call MPI_Init(ierror)
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierror)
   dims = [0, 0]
   call MPI_Dims_Create(np, 2_int32, dims)
   npcol = np / dims(2)
   nprow = dims(2)

   if (myrank == 0) then 
      write(stdout,*) 'Proc = ', np
      write(stdout,*) 'Rows = ', nprow
      write(stdout,*) 'Cols = ', npcol
   endif

   call blacs_get(0, 0, ictxt)
   call blacs_gridinit(ictxt, 'row', nprow , npcol)
   call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
   call initialize_matrix(C, N1 + N2, N1 + N2, Cdesc, init_val = ONE)
   call initialize_matrix(B, N2     , N2     , Bdesc, init_val = N1*ONE)

   alpha = -ONE 
   beta  =  ONE 
  
   call PGEMM ('N', 'N', N2, N2, N1, &
      alpha, C, N1 + 1, 1, Cdesc, C, 1, N1 + 1, Cdesc, &
      beta, B, 1, 1, Bdesc)

   call check_matrix(B, Bdesc)

   !call print_matrix(B, N2, N2, 1, 1, Bdesc, '/tmp/B')
   
   call blacs_gridexit(ictxt)
   call blacs_exit(1) 
   call MPI_Finalize(ierror)

contains

   subroutine initialize_matrix(A, M, N, desc, init_val)
#if defined(fD) 
      real(kind=dp), allocatable, intent(inout) :: A(:)
      real(kind=dp), intent(in), optional :: init_val
#elif defined(fC)   
      complex(kind=sp), allocatable, intent(inout) :: A(:)
      complex(kind=sp), intent(in), optional :: init_val
#elif defined(fZ)   
      complex(kind=dp), allocatable, intent(inout) :: A(:)
      complex(kind=dp), intent(in), optional :: init_val
#else
      real(kind=sp), allocatable, intent(inout) :: A(:)
      real(kind=sp), intent(in), optional :: init_val
#endif
      integer, intent(in) :: M, N
      integer, intent(out) :: desc(dlen_)

      integer :: mloc, nloc, lld, info

      desc(ctxt_) = -1
      if (allocated(A)) deallocate(A)

      if (myprow >= 0) then
         mloc = numroc(M, mb, myprow, 0, nprow)
         nloc = numroc(N, nb, mypcol, 0, npcol)
         lld = max(1, mloc)
         call descinit(desc, M, N, mb, nb, 0, 0, ictxt, lld, info)
         allocate(A(lld*nloc))
         if (present(init_val)) A = init_val
      endif
   end subroutine initialize_matrix
   
   subroutine check_matrix(A, desc)
#if defined(fD) 
      real(kind=dp), intent(in) :: A(:)
#elif defined(fC) 
      complex(kind=sp), intent(in) :: A(:)
#elif defined(fZ) 
      complex(kind=dp), intent(in) :: A(:)
#else
      real(kind=sp), intent(in) :: A(:)
#endif
      integer, intent(in) :: desc(dlen_)

      integer :: ir, ic, mloc, nloc, counter, idx

      counter = 0
      
      if (myprow >= 0) then
         mloc = numroc(desc(m_), desc(mb_), myprow, 0, nprow)
         nloc = numroc(desc(n_), desc(nb_), mypcol, 0, npcol)

         do ir = 1, mloc
            do ic = 1, nloc
               idx = (ic - 1)*desc(lld_) + ir
               if (A(idx) /= ZERO) then 
                  print *, &
                     indxl2g(ir , mb, myprow, 0, nprow), &
                     indxl2g(ic , nb, mypcol, 0, npcol), A(idx)
                  counter = counter + 1
               endif
            enddo
         enddo
      endif
   end subroutine check_matrix

end program test_pzgemm
