
!==========================================================
!==========================================================
!==========================================================
program proj_qr_svd_scalapack
#if defined(_OPENMP)
  use omp_lib
#endif
  use low_level
  implicit none

  real(dp), allocatable :: A(:,:), Y(:,:), B(:,:), C(:,:)
  real(dp), allocatable :: tau(:)
  integer :: descA(NDEL), descY(NDEL), descB(NDEL), descC(NDEL)

  integer :: kp
  integer :: ifile
  character(len=128) :: input_file_name
  !=====

#if defined(_OPENMP)
  write(*,*) 'OPENMP threads:', omp_get_max_threads()
#endif

  call init_scalapack()

  ! Default
  nI = 73**2
  nG = 400 * 2
  method = 'PROJ_QR_SVD'
  k = 400
  p = 0
  q = 1

  if( COMMAND_ARGUMENT_COUNT() == 1 ) then
    call GET_COMMAND_ARGUMENT(1, VALUE=input_file_name)
  else
    stop "no input file"
  endif

  if( rank == 0 ) write(stdout,*) 'Reading input file: ', input_file_name
  open(newunit=ifile, file=TRIM(input_file_name), status='old', action='read')
  read(ifile,input)
  close(ifile)

  ! Assume nI >> nG
  ! Enforce it
  if( nI <= nG ) stop "nI <= nG"

  if( rank == 0 ) write(*,*) 'k=',k
  if( rank == 0 ) write(*,*) 'p=',p
  if( rank == 0 ) write(*,*) 'q=',q

  kp = k + p

  call get_matrix_A(file_in, nI, nG, A, descA)

  if( rank == 0 ) write(*,*) 'Proc', rank,'local A', SIZE(A,DIM=1), SIZE(A,DIM=2) 
  if( rank == 0 ) write(*,*) 'Proc', rank,'global A', descA(M_), descA(N_)


  select case(TRIM(method))
  case('PROJ_QR_SVD')
    ! proj
    call step1(kp, q, A, descA, Y, descY)
    ! QR
    call step2(Y, descY, tau)
    ! SVD
    call step3(Y, descY, tau, A, descA, B, descB)
    call step4(Y, descY, B, descB, C, descC)
    call step5(Y, descY, tau, C, descC)
  case('SVD')
    call full_svd(A, descA, C, descC)
  case default
    stop "method not valid"
  end select

  !!DEBUG
  !block
  !  real(dp) :: rtmp
  !  call PDDOT( k, rtmp, C, 1, 1, descC, nI, C, 1, 1, descC, nI)
  !  write(stdout,*) 'Proc', rank, 'final integral (11|11) (Ha):', rtmp
  !end block

  call dump_matrix_C(k, file_out, C, descC)
  call finalize_scalapack()

  if( rank == 0 ) write(*,*) 'Job done'


end program proj_qr_svd_scalapack


