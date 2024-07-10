!==========================================================
program proj_qr_svd_scalapack
#if defined(_OPENMP)
  use omp_lib
#endif
  use low_level
  use input_param
  implicit none

  real(dp), allocatable :: A(:,:), Y(:,:), B(:,:), C(:,:)
  real(dp), allocatable :: tau(:)
  integer :: descA(NDEL), descY(NDEL), descB(NDEL), descC(NDEL)

  !=====


  call init_scalapack()

#if defined(_OPENMP)
  if( rank == 0 ) write(stdout,*) 'OPENMP threads:', omp_get_max_threads()
#endif


  call read_input_file()


  call get_matrix_A(file_in, nI, nG, A, descA)

  if( rank == 0 ) write(stdout,*) 'Proc', rank,'local A', SIZE(A,DIM=1), SIZE(A,DIM=2) 
  if( rank == 0 ) write(stdout,*) 'Proc', rank,'global A', descA(M_), descA(N_)


  select case(TRIM(method))
  case('PROJ_QR_SVD')
    ! Random noise to guess Span(A)
    call step1(kp, q, A, descA, Y, descY)
    ! QR decomposition of Y
    call step2(Y, descY, tau)
    ! Apply Q on A -> B
    call step3(Y, descY, tau, A, descA, B, descB)
    ! SVD of B
    call step4(Y, descY, B, descB, C, descC)
    ! Apply Q on C
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

  if( rank == 0 ) write(stdout,*) 'Job done'


end program proj_qr_svd_scalapack


