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
   integer, parameter :: MAX_PRINTED_LINES_PER_PROC=5

#if defined(fD)
#define PGEMM pdgemm
#define PGERU pdger
#define PSYRK pdsyrk
   real(kind=dp), allocatable :: C(:), B(:)
   real(kind=dp) :: alpha, beta
   real(kind=dp), parameter :: ONE = 1.0_dp, ZERO = 0.0_dp
#elif defined(fC) 
#define PGEMM pcgemm
#define PGERU pcgeru
#define PSYRK pcsyrk
   complex(kind=sp), allocatable :: C(:), B(:)
   complex(kind=sp) :: alpha, beta
   complex(kind=sp), parameter :: ONE = (1.0_sp, 0.0_sp), ZERO = (0.0_sp, 0.0_sp)
#elif defined(fZ)   
#define PGEMM pzgemm
#define PGERU pzgeru
#define PSYRK pzsyrk
   complex(kind=dp), allocatable :: C(:), B(:)
   complex(kind=dp) :: alpha, beta
   complex(kind=dp), parameter :: ONE = (1.0_dp, 0.0_dp), ZERO = (0.0_dp, 0.0_dp)
#else
#define PGEMM psgemm
#define PGERU psger
#define PSYRK pssyrk
   real(kind=sp), allocatable :: C(:), B(:)
   real(kind=sp) :: alpha, beta
   real(kind=sp), parameter :: ONE = 1.0_sp, ZERO = 0.0_sp
#endif

   integer :: Cdesc(dlen_), Bdesc(dlen_)
   integer :: ictxt, nprow, npcol, myprow, mypcol

   integer(kind=int32) :: np, ierror, myrank, dims(2)
   integer :: i, nfails, total_nfails

   call MPI_Init(ierror)
   call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierror)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierror)
   dims = [0, 0]
   call MPI_Dims_Create(np, 2_int32, dims)
   npcol = np / dims(2)
   nprow = dims(2)

   if (myrank == 0) then 
      write(stdout, '(/,"Procs = ", I0, "; Grid = (", I0, " x ", I0, ")")') &
         np, nprow, npcol
      write(stdout, '("mb = ", I0)', advance='no') mb
      write(stdout, '("; nb = ", I0)', advance='no') nb
      write(stdout, '("; N1 = ", I0)', advance='no') N1
      write(stdout, '("; N2 = ", I0)') N2
   endif

   call blacs_get(0, 0, ictxt)
   call blacs_gridinit(ictxt, 'row', nprow , npcol)
   call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
   call initialize_matrix(C, N1 + N2, N1 + N2, Cdesc, init_val = ONE)
   call initialize_matrix(B, N2     , N2     , Bdesc, init_val = N1*ONE)

   alpha = -ONE 
   beta  =  ONE 
  
   ! Test using p?gemm
   call PGEMM ('N', 'N', N2, N2, N1, &
      alpha, C, N1 + 1, 1, Cdesc, C, 1, N1 + 1, Cdesc, &
      beta, B, 1, 1, Bdesc)
   ! Test using p?ger or p?geru
   !do i = 1, N1
   !   call PGERU (N2, N2, &
   !      alpha, &
   !      C, N1 + 1, i, Cdesc, 1, &
   !      C, i, N1 + 1, Cdesc, N1+N2, &
   !      B, 1, 1, Bdesc)
   !enddo

   ! Test using p?syrk
   !call PSYRK ('L', 'N', N2, N1, &
   !   alpha, C, N1 + 1, 1, Cdesc, &
   !   beta, B, 1, 1, Bdesc)

   call check_matrix(B, Bdesc, nfails)
   call MPI_Reduce(nfails, total_nfails, 1_int32, MPI_INTEGER, MPI_SUM, 0_int32, MPI_COMM_WORLD, ierror)
   if (myrank == 0) then
      if (total_nfails == 0) write (stdout, '(A)') "TEST PASSED!"
   endif

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
   
   subroutine check_matrix(A, desc, counter)
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
      integer, intent(out) :: counter

      integer :: lr, lc, gr, gc, mloc, nloc, idx

      counter = 0
      
      if (myprow >= 0) then
         mloc = numroc(desc(m_), desc(mb_), myprow, 0, nprow)
         nloc = numroc(desc(n_), desc(nb_), mypcol, 0, npcol)

         do lr = 1, mloc
            do lc = 1, nloc
               idx = (lc - 1)*desc(lld_) + lr
               gr = indxl2g(lr , mb, myprow, 0, nprow)
               gc = indxl2g(lc , nb, mypcol, 0, npcol)
               !if (gr < gc) cycle
               if (A(idx) /= ZERO) then 
                  if (counter < MAX_PRINTED_LINES_PER_PROC) then
                     write (stdout, '("B(", I0, ", ", I0, ")", T15, " = ")', advance="no") gr, gc
                     write (stdout, *) A(idx)
                  endif
                  counter = counter + 1
               endif
            enddo 
         enddo
         if (counter >= MAX_PRINTED_LINES_PER_PROC) write(stdout, '(A)') '...'
      endif
      call blacs_barrier(ictxt, 'A')
      if (counter > 0) write(stdout, &
            '("TEST FAILED in proc (", I0, ", ", I0, ") with ", I0, " errors")') &
            myprow, mypcol, counter
   end subroutine check_matrix

end program test_pzgemm
