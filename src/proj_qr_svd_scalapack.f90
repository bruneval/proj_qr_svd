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

#if defined(_OPENMP)
  write(*,*) 'OPENMP threads:', omp_get_max_threads()
#endif

  call init_scalapack()


  call read_input_file()


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


