program qr_svd
  use input_param
  implicit none

  !integer, parameter :: stdout=OUTPUT_UNIT
  !integer, parameter :: dp=8
  integer :: m, n
  real(dp), allocatable :: A(:,:), B(:,:), C(:,:)
  real(dp), allocatable ::  M_ref(:,:), M_qrsvd(:,:), M_svd(:,:), M_projqrsvd(:,:)
  real(dp), allocatable :: R(:,:)
  real(dp), allocatable :: U(:,:)
  real(dp), allocatable :: VT(:,:)
  real(dp), allocatable :: Omega(:,:), Y(:,:)
  real(dp), allocatable :: tau(:), sigma(:)
  real(dp),allocatable :: work(:)
  integer :: lwork, info, i, j, iq, ifile
  real(dp) :: start, finish, start0, finish0
  integer :: descA(NDEL), descC(NDEL)
  !character(len=128) :: method, file_out, file_in

  !file_in  = "lih_fc_fv_irene_CoulombVertex.elements"
  !file_out = "test.elements"
  !method = "QR_SVD"
  !method = "SVD"
  !method = "PROJ_QR_SVD"
  !k = 500
  !q = 1

  call init_scalapack()
  call read_input_file()

  write(*,*) 'k=',k
  write(*,*) 'q=',q
  write(*,*) 'nproc=',nproc
  write(*,*) nG, npw, nI

  call get_matrix_A(file_in, nI, nG, A, descA)
  m = SIZE(A,dim=1)
  n = SIZE(A,dim=2)
  write(*,*) 'sizes m, n', m ,n
  write(*,*) 'A in Mb', 8.0 * m * n / 1024.**2
  allocate( VT(1,1) )

!  write(*,*) 
!  write(*,*)  '================= REFERENCE ========================='
!
!  allocate(M_ref(m,m))
!  call cpu_time(start)
!  call DGEMM('N', 'T', m, m, n, 1.0d0, A, m, A, m, 0.0d0, M_ref, m)
!  !call DSYRK('U', 'N', m, n, 1.0d0, A, m, 0.0d0, M_ref, m)
!  call cpu_time(finish)
!  write(*,*) 'DGEMM time:', finish - start, 'seconds'
!  write(*,*) '(1|1)', M_ref(1,1)
!  write(*,*) '(2|2)', M_ref(2,2)
!  write(*,*) '(1|2)', M_ref(1,2)
!  write(*,*) '(2|1)', M_ref(2,1)
!  deallocate(M_ref)

 
  
  ! Assume m >> n
  if( m <= n ) stop
  

select case(TRIM(method))
case('SVD')
  !
  ! Full SVD
  write(*,*) 
  write(*,*)  '================= FULL SVD ========================='
  call cpu_time(start)
  allocate( sigma(n) )
  allocate( U(m,n) )
  allocate(work(1))
  lwork = -1
  U(:,:) = 0.0d0
  call DGESVD("S", "N", m, n, A, m, sigma, U, m, VT, 1, work, lwork, info)
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  call DGESVD("S", "N", m, n, A, m, sigma, U, m, VT, 1, work, lwork, info)
  deallocate(work)

  allocate( C(m,k) )
  do i=1,k
    C(:,i) = U(:,i) * sigma(i)
    write(100,*) sigma(i)
  enddo
  call cpu_time(finish)
  write(*,*) 'full SVD time:', finish - start, 'seconds'

!  allocate(M_svd(m,m))
!  call cpu_time(start)
!  call DGEMM('N', 'T', m, m, k, 1.0d0, U, m, U, m, 0.0d0, M_svd, m)
!  call cpu_time(finish)
!  write(*,*) 'DGEMM time:', finish - start, 'seconds'
!  write(*,*) '(1|1)', M_svd(1,1)
!  write(*,*) '(2|2)', M_svd(2,2)
!  write(*,*) '(1|2)', M_svd(1,2)
!  write(*,*) '(2|1)', M_svd(2,1)
!  deallocate(M_svd)


case('QR_SVD')
  !
  ! QR SVD
  write(*,*) 
  write(*,*)  '================= QR SVD ========================='

  call cpu_time(start)
  allocate(tau(n))
  allocate(R(n,n))
  allocate(work(1))
  lwork = -1
  call DGEQRF(m, n, A, m, tau, work, lwork, info)
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  call DGEQRF(m, n, A, m, tau, work, lwork, info)
  deallocate(work)

  R(:,:) = 0.0_dp
  do i = 1, m
    do j = 1, n
      if (i <= j) then
        R(i,j) = A(i,j)
      endif
    enddo
  enddo
  call cpu_time(finish)
  write(*,*) 'QR time:', finish - start, 'seconds'


  call cpu_time(start)
  allocate(work(1))
  lwork = -1
  U(:,:) = 0.0d0
  call DGESVD("S", "N", n, n, R, n, sigma, U, m, VT, 1, work, lwork, info)
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  call DGESVD("S", "N", n, n, R, n, sigma, U, m, VT, 1, work, lwork, info)
  deallocate(work)

  do i=1,n
    U(:,i) = U(:,i) * sigma(i)
  enddo
  deallocate(sigma)

  call cpu_time(finish)
  write(*,*) 'small SVD time:', finish - start, 'seconds'

  ! DORMQR applies Q on a matrix U
  call cpu_time(start)
  allocate(work(1))
  lwork = -1
  call DORMQR( "L", "N", m, k, n, A, m, tau, U, m, work, lwork, info)
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  call DORMQR( "L", "N", m, k, n, A, m, tau, U, m, work, lwork, info)
  deallocate(work)
  deallocate(tau)
  deallocate(R)

  call cpu_time(finish)
  write(*,*) 'Q*U product time:', finish - start, 'seconds'

  !allocate(M_qrsvd(m,m))
  !call cpu_time(start)
  !call DGEMM('N', 'T', m, m, k, 1.0d0, U, m, U, m, 0.0d0, M_qrsvd, m)
  !call cpu_time(finish)
  !deallocate(U)
  !write(*,*) 'DGEMM time:', finish - start, 'seconds'
  !write(*,*) '(1|1)', M_qrsvd(1,1)
  !write(*,*) '(2|2)', M_qrsvd(2,2)
  !write(*,*) '(1|2)', M_qrsvd(1,2)
  !write(*,*) '(2|1)', M_qrsvd(2,1)
  !deallocate(M_qrsvd)

  allocate(C(m,k))
  C(:,:) = U(:,1:k)

case('PROJ_QR_SVD')

  !
  ! proj QR SVD
  write(*,*) 
  write(*,*)  '================= proj QR SVD ========================='

  
  call cpu_time(start0)

  allocate(Omega(n,k))
  allocate(Y(m,k))

  !
  ! Step 1: Create Y
  !
  write(*,*) ' **** Step 1 **** '
  call cpu_time(start)

  ! Random Omega centered in zero
  if( .TRUE. ) then
    call random_number(Omega)
    Omega(:,:) = Omega(:,:) - 0.50d0
  else
    write(*,*) 'Reading Omega from Omega.dat'
    open(newunit=ifile,file='Omega.dat',action='read')
    do i = 1, n
      do j = 1, k
        read(ifile,*) Omega(i,j)
      enddo
    enddo
    close(ifile)
  endif

  ! Y = A * Omega
  call DGEMM('N', 'N', m, k, n, 1.0d0, A, m, Omega, n, 0.0d0, Y, m)
  do iq=1,q
    ! Omega = A**T * Y
    call DGEMM('T', 'N', n, k, m, 1.0d0, A, m, Y, m, 0.0d0, Omega, n)
    ! Y = A * Omega
    call DGEMM('N', 'N', m, k, n, 1.0d0, A, m, Omega, n, 0.0d0, Y, m)
  enddo
  call cpu_time(finish)
  write(*,*) 'DGEMMs from noisy Omega time:', finish - start, 'seconds'

  !
  ! Step 2: Q R decomposition of Y
  !
  write(*,*) ' **** Step 2 **** '
  call cpu_time(start)
  allocate(tau(k))
  allocate(work(1))
  lwork = -1
  call DGEQRF(m, k, Y, m, tau, work, lwork, info)
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  call DGEQRF(m, k, Y, m, tau, work, lwork, info)
  deallocate(work)

  ! Get R from upper triangle of output Y
  allocate(R(k,k))
  R(:,:) = 0.0_dp
  do i = 1, m
    do j = 1, k
      if (i <= j) then
        R(i,j) = Y(i,j)
      endif
    enddo
  enddo
  call cpu_time(finish)
  write(*,*) 'QR time:', finish - start, 'seconds'


  !
  ! Step 3: Create B
  !
  ! B = Q**T * A
  allocate( B(k,n) )
  write(*,*) ' **** Step 3 **** '

  ! DORMQR applies Q**T on a matrix A
  call cpu_time(start)
  allocate(work(1))
  lwork = -1
  call DORMQR( "L", "T", m, n, k, Y, m, tau, A, m, work, lwork, info)
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  call DORMQR( "L", "T", m, n, k, Y, m, tau, A, m, work, lwork, info)
  deallocate(work)
  B(1:k,:) = A(1:k,:)

  call cpu_time(finish)
  write(*,*) 'B = Q**T*A product time:', finish - start, 'seconds'
  deallocate(R)

  !call print_matrix('B_fortan_lapack.dat',A)

  write(*,*) ' ************* B ************* '
  write(*,'(*(f12.4,1x))') B(1,1:5)
  write(*,'(*(f12.4,1x))') B(2,1:5)
  write(*,'(*(f12.4,1x))') B(3,1:5)
  write(*,'(*(f12.4,1x))') B(4,1:5)
  write(*,'(*(f12.4,1x))') B(5,1:5)


  !
  ! Step 4: SVD of B
  !
  write(*,*) ' **** Step 4 **** '

  allocate(sigma(k))
  allocate(U(k,k))

  call cpu_time(start)
  allocate(work(1))
  lwork = -1
  !U(:,:) = 0.0d0
  call DGESVD("S", "N", k, n, B, k, sigma, U, k, VT, 1, work, lwork, info)
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  call DGESVD("S", "N", k, n, B, k, sigma, U, k, VT, 1, work, lwork, info)
  deallocate(work)

  do i=1,k
    U(:,i) = U(:,i) * sigma(i)
    write(101,*) sigma(i)
  enddo
  call cpu_time(finish)
  deallocate(sigma)
  write(*,*) 'B SVD time:', finish - start, 'seconds'

  ! Step 5: Apply Q on U
  write(*,*) ' **** Step 5 **** '
  ! DORMQR applies Q on a matrix U
  call cpu_time(start)
  allocate(C(m,k))
  C(:,:) = 0.0d0
  C(1:k,:) = U(1:k,:)
  allocate(work(1))
  lwork = -1
  call DORMQR( "L", "N", m, k, k, Y, m, tau, C, m, work, lwork, info)
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  call DORMQR( "L", "N", m, k, k, Y, m, tau, C, m, work, lwork, info)
  deallocate(work)
  deallocate(tau)

  call cpu_time(finish)
  write(*,*) 'C = Q*U product time:', finish - start, 'seconds'

  call cpu_time(finish0)
  write(*,*) 'Proj QR SVD total time:', finish0 - start0, 'seconds'

  !! Final check
  !allocate(M_projqrsvd(m,m))
  !call cpu_time(start)
  !call DGEMM('N', 'T', m, m, k, 1.0d0, C, m, C, m, 0.0d0, M_projqrsvd, m)
  !call cpu_time(finish)
  !deallocate(U)
  !write(*,*) 'DGEMM time:', finish - start, 'seconds'
  !write(*,*) '(1|1)', M_projqrsvd(1,1)
  !write(*,*) '(2|2)', M_projqrsvd(2,2)
  !write(*,*) '(1|2)', M_projqrsvd(1,2)
  !write(*,*) '(2|1)', M_projqrsvd(2,1)
  !deallocate(M_projqrsvd)


case default
  stop "error in method choice"
end select

  ! Fake descriptor
  call DESCINIT(descC, nI, kp, block_row, block_col, first_row, first_col, cntxt_sd, nI, info)

  call dump_matrix_C(k,file_out, C, descC)



contains
subroutine getget_matrix_A(file_in, A)
    implicit none
    character(len=*),intent(in) :: file_in
    real(dp),allocatable,intent(inout) :: A(:,:)
    logical,parameter :: random_matrix = .FALSE.
    integer :: I, G
    integer :: unitcv, nstate, istate, jstate, ijstate, npw
    integer :: nstate_file, ijstate_file
    complex(dp),allocatable :: coulomb_vertex_ij(:)
 
    if( random_matrix ) then
        write(*,'(/,1x,a)') 'Random matrix for testing'
        I = 1000   ! Number of rows (can be millions in actual scenario)
        G = 800    ! Number of columns (can be large)
        allocate(A(I, G))
        ! Initialize matrix A with random values
        call random_number(A)
    else
        write(*,*) 'Opening file:', TRIM(file_in)
        ! hard-coded
        nstate_file= 172
        nstate = 172
        npw = 6662
        allocate(coulomb_vertex_ij(npw))
        I = nstate**2
        G = npw * 2
        allocate(A(I, G))

        open(newunit=unitcv, file=TRIM(file_in), form='unformatted', access='stream', status='old', action='read')
        ijstate = 0
        do istate=1,nstate_file
          do jstate=1,nstate_file
            read(unitcv) coulomb_vertex_ij(:)
            if( istate > nstate .OR. jstate > nstate ) cycle
            ijstate = ijstate + 1
            
            A(ijstate,1:npw)       = coulomb_vertex_ij(:)%re
            A(ijstate,npw+1:2*npw) = coulomb_vertex_ij(:)%im
          enddo
        enddo
      
        close(unitcv)
        write(*,*) 'Closing file:', TRIM(file_in)

    endif

end subroutine getget_matrix_A


subroutine dumpdump_matrix_C(file_out, C)
  implicit none

  character(len=*),intent(in) :: file_out
  real(dp),allocatable,intent(inout)  :: C(:,:)
  !=====
  integer :: nI, k
  integer :: unitcv, ierr, info
  complex(dp),allocatable :: coulomb_vertex_I(:)
  integer :: kc, Ig
  !=====

  nI = SIZE(C,DIM=1)
  k  = SIZE(C,DIM=2)
  kc = k / 2


  write(stdout,'(/,1x,a,a,a)') 'Writing file ', TRIM(file_out), ' with fortran'
  write(stdout,'(5x,a,i6,a,i8)') 'CoulombVertex dimensions:', kc, ' x ', nI

  allocate(coulomb_vertex_I(kc))
  write(*,*) STORAGE_SIZE(coulomb_vertex_I(1)) / 8.0
  write(*,*) STORAGE_SIZE(coulomb_vertex_I(1)) * kc * nI / 1024.**2 /8.0

  open(newunit=unitcv, file=TRIM(file_out), access='stream', action='write')

  do Ig=1, nI
    coulomb_vertex_I(:) = CMPLX(C(Ig,1:kc), C(Ig,kc+1:2*kc))

    write(unitcv) coulomb_vertex_I(:)

  enddo
  close(unitcv)

end subroutine dumpdump_matrix_C

subroutine print_matrix(name, matrix)
  implicit none
  character(len=*),intent(in) :: name
  real(dp) :: matrix(:,:)
  integer :: m, n, i, j, ifile

  m = SIZE(matrix,DIM=1)
  n = SIZE(matrix,DIM=2)
  open(newunit=ifile,file=TRIM(name),action='write')
  do i=1,m
    do j=1,n
      write(ifile,'(es12.3)') matrix(i,j)
    enddo
  enddo
  close(ifile)
end subroutine print_matrix

end program
