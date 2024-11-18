! 
! The program test pdgemm matrix x matrix multiplication under fixed condition 
! on a square processor grid provided by the user.
! 
! The product tested is:
! 
!  C = A * B
! 
! with A being a 8160 x 8160  matrix with all coeffs set to 1
! and  B being a 8160 x 19140 matrix with all coeffs set to 1
! The result expected is thus all coeffs of C equal to 8160
! 
PROGRAM TEST
  
  ! Parameters
  INTEGER         , PARAMETER :: M=8160, N =19140, K=8160, DLEN_=9
  INTEGER         , PARAMETER :: CSRC=1, RSRC=1
  DOUBLE PRECISION, PARAMETER :: ONE=1.0D+0, ZERO=0.0D+0
  
  ! work variables
  INTEGER                                       :: ICTXT
  INTEGER                                       :: IAM
  INTEGER                                       :: NPROCS
  INTEGER                                       :: NPROW
  INTEGER                                       :: NPCOL
  INTEGER                                       :: MYROW
  INTEGER                                       :: MYCOL
  INTEGER                                       :: DESCA(9)
  INTEGER                                       :: DESCB(9)
  INTEGER                                       :: DESCC(9)
  INTEGER                                       :: M_A
  INTEGER                                       :: N_A
  INTEGER                                       :: M_B
  INTEGER                                       :: N_B
  INTEGER                                       :: M_C
  INTEGER                                       :: N_C
  INTEGER                                       :: MB_A
  INTEGER                                       :: NB_A
  INTEGER                                       :: MB_B
  INTEGER                                       :: NB_B
  INTEGER                                       :: MB_C
  INTEGER                                       :: NB_C
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: B
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: C
  
  ! Get starting information
  CALL BLACS_PINFO( IAM, NPROCS )
  
  ! try setting square grid
  NPROW = sqrt(REAL(NPROCS,kind=8))
  NPCOL = sqrt(REAL(NPROCS,kind=8))
  if ( NPROW*NPCOL .ne. NPROCS ) then
    print *,"please provide a square number of procs"
    stop 1
  end if
  
  ! Define process grid
  CALL BLACS_GET( -1, 0, ICTXT )
  CALL BLACS_GRIDINIT( ICTXT, 'R', NPROW, NPCOL )
  CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
  
  ! set A matrix dimensions
  M_A = M
  N_A = K
  
  ! set B matrix dimensions
  M_B = K
  N_B = N
  
  ! set C matrix dimensions
  M_C = M
  N_C = N
  
  ! set blocking factors for A matrix
  MB_A = M_A/NPROW
  NB_A = N_A/NPCOL
  
  ! set blocking factors for B matrix
  MB_B = M_B/NPROW
  NB_B = 32
  
  ! set blocking factors for C matrix
  MB_C = M_C/NPROW
  NB_C = 32
  
  ! get A local dimensions
  MLOC_A = NUMROC( M_A, MB_A, MYROW, 0, NPROW )
  NLOC_A = NUMROC( N_A, NB_A, MYCOL, 0, NPCOL )
  
  ! get B local dimensions
  MLOC_B = NUMROC( M_B, MB_B, MYROW, 0, NPROW )
  NLOC_B = NUMROC( N_B, NB_B, MYCOL, 0, NPCOL )
  
  ! get C local dimensions
  MLOC_C = NUMROC( M_C, MB_C, MYROW, 0, NPROW )
  NLOC_C = NUMROC( N_C, NB_C, MYCOL, 0, NPCOL )
  
  ! Initialize the array descriptor for the matrix A, B and C
  CALL DESCINIT( DESCA, M_A, N_A, MB_A, NB_A, 0, 0, ICTXT, max(MLOC_A,1), INFO )
  CALL DESCINIT( DESCB, M_B, N_B, MB_B, NB_B, 0, 0, ICTXT, max(MLOC_B,1), INFO )
  CALL DESCINIT( DESCC, M_C, N_C, MB_C, NB_C, 0, 0, ICTXT, max(MLOC_C,1), INFO )
  
  ! print grid infos
  do IPROC=0,NPROCS-1
    if ( IPROC .eq. IAM ) then
      print *,""
      print *,"-------------------------"
      print *,"PROC, MYROW, MYCOL :",PROC,MYROW,MYCOL
      print *,"MLOC_A, NLOC_A :",MLOC_A,NLOC_A
      print *,"MLOC_B, NLOC_B :",MLOC_B,NLOC_B
      print *,"MLOC_C, NLOC_C :",MLOC_C,NLOC_C
      print *,"DESCA :",DESCA
      print *,"DESCB :",DESCB
      print *,"DESCC :",DESCC
      print *,"-------------------------"
      print *,""
    end if
    CALL SLEEP(2)
  end do
  
  ! allocate and set matrices
  ALLOCATE( A(MLOC_A,NLOC_A) )
  ALLOCATE( B(MLOC_B,NLOC_B) )
  ALLOCATE( C(MLOC_C,NLOC_C) )
  
  ! init A matrix
  do j=1,NLOC_A
    do i=1,MLOC_A
      A(i,j)=ONE
    end do
  end do
  
  ! init B matrix
  do j=1,NLOC_B
    do i=1,MLOC_B
      B(i,j)=ONE
    end do
  end do
  
  ! compute A * B
  CALL PDGEMM('N', 'N',       &
&             M, N, K,        &
&             ONE,            &
&             A, 1, 1, DESCA, &
&             B, 1, 1, DESCB, &
&             ZERO,           &
&             C, 1, 1, DESCC )
  
  ! check result
  do j=1,NLOC_C
    do i=1,MLOC_C
      if ( abs(C(i,j)-K) .gt. 1.0D-8 ) then
        print *,"Error: result differs from exact"
        print *,"C(",i,",",j,")=",C(i,j)
        print *,"expected ",K
        print *,"TEST FAILED!"
        stop 2
      end if
    end do
  end do
 
  ! inform that everything is ok
  if ( IAM .eq. 0 ) then
    print *,"TEST PASSED!"
  end if 

  ! terminate
  CALL BLACS_GRIDEXIT( ICTXT )
  CALL BLACS_EXIT( 0 )
  
END PROGRAM
