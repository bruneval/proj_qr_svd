module low_level
  use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use mpi

  integer,parameter :: stdout = OUTPUT_UNIT
  integer,parameter :: dp = 8
  integer,parameter :: NDEL = 9
  integer,parameter :: M_ = 3
  integer,parameter :: N_ = 4
  integer,parameter :: SCALAPACK_BLOCKSIZE_MAX = 64
  integer,parameter :: block_row = SCALAPACK_BLOCKSIZE_MAX
  integer,parameter :: block_col = SCALAPACK_BLOCKSIZE_MAX
  integer,parameter :: first_row = 0
  integer,parameter :: first_col = 0
  integer,external  :: NUMROC, INDXL2G, INDXG2L, INDXG2P

  integer :: rank, nproc
  integer :: cntxt_cd, nprow_cd, npcol_cd, iprow_cd, ipcol_cd
  integer :: cntxt_sd, nprow_sd, npcol_sd, iprow_sd, ipcol_sd
  integer :: cntxt_one, nprow_one, npcol_one, iprow_one, ipcol_one


contains


subroutine init_scalapack()
  implicit none

  integer :: info
  !=====

  call MPI_INIT(info)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, info)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, info)

  ! Squared grid
  nprow_sd = INT(SQRT(REAL(nproc)))
  npcol_sd = nproc / nprow_sd
  do while( npcol_sd <= nproc .AND. nprow_sd >= 1 )
    ! Found a correct distribution
    if( nprow_sd * npcol_sd == nproc ) exit
    npcol_sd = npcol_sd + 1
    nprow_sd = nproc / npcol_sd
  enddo

  !!!DEBUG
  !nprow_sd = nproc
  !npcol_sd = 1

  if( nprow_sd * npcol_sd /= nproc ) stop "Error in the squared grid"
  call BLACS_GET( -1, 0, cntxt_sd )
  call BLACS_GRIDINIT( cntxt_sd, 'R', nprow_sd, npcol_sd )
  call BLACS_GRIDINFO( cntxt_sd, nprow_sd, npcol_sd, iprow_sd, ipcol_sd )

  ! Column-only grid
  nprow_cd = 1
  npcol_cd = nproc
  call BLACS_GET( -1, 0, cntxt_cd )
  call BLACS_GRIDINIT( cntxt_cd, 'R', nprow_cd, npcol_cd )
  call BLACS_GRIDINFO( cntxt_cd, nprow_cd, npcol_cd, iprow_cd, ipcol_cd )

  ! Fake cntxt with one proc
  nprow_one= 1
  npcol_one= 1
  call BLACS_GET( -1, 0, cntxt_one)
  call BLACS_GRIDINIT( cntxt_one, 'R', nprow_one, npcol_one)
  call BLACS_GRIDINFO( cntxt_one, nprow_one, npcol_one, iprow_one, ipcol_one)

end subroutine init_scalapack


subroutine finalize_scalapack()
  implicit none

  integer :: ier
  !=====

  call BLACS_GRIDEXIT( cntxt_cd )
  call BLACS_GRIDEXIT( cntxt_sd )
  call MPI_FINALIZE(ier)

end subroutine finalize_scalapack


subroutine get_matrix_A(file_in, nI, nG, At, descAt)
  implicit none

  character(len=*),intent(in) :: file_in
  integer,intent(in) :: nI, nG
  real(dp),allocatable,intent(out) :: At(:,:)
  integer,intent(out) :: descAt(NDEL)
  !=====
  real(dp),allocatable :: Aread(:,:), A(:,:)
  integer :: npw
  integer :: unitcv, ierr
  complex(dp),allocatable :: coulomb_vertex_I(:)
  integer :: complex_length
  real(dp) :: rtmp
  integer(kind=MPI_OFFSET_KIND) :: disp, disp_increment
  integer :: descAread(NDEL), descA(NDEL)
  real(dp), allocatable :: eri_3center_tmp(:,:)
  integer :: mAr, nAr, mA, nA, mAt, nAt, info
  integer :: nstate2
  integer :: Ig, Il

  write(stdout,'(1x,a,i7,a,i7)') 'Dimensions read:', nG, ' x ', nI

  npw = nG / 2
  if( npw * 2 /= nG ) stop 'nG should be even'
  allocate(coulomb_vertex_I(npw))
  write(*,*) 'npw',npw,nG

  ! Create a SCALAPACK matrix (nG, nI) that is distributed on column index only
  mAr = NUMROC(nG, block_row, iprow_cd, first_row, nprow_cd)
  nAr = NUMROC(nI, block_col, ipcol_cd, first_col, npcol_cd)
  
  if( rank == 0 ) write(stdout,'(/,1x,a,a,a)') 'Reading file ', TRIM(file_in), ' with MPI-IO'
  if( rank == 0 ) write(stdout,'(5x,a30,i4,a,i4)') 'using a processor grid:', nprow_cd, ' x ', npcol_cd
  if( rank == 0 ) write(stdout,'(5x,a30,i6,a,i8)') 'CoulombVertex dimensions:', npw, ' x ', nI

  if( nproc > 1 ) then

    call DESCINIT(descAread, nG, nI, block_row, block_col, first_row, first_col, cntxt_cd, mAr, info)
    if( info /= 0 ) stop "DESCINIT failure"

    allocate(Aread(mAr,nAr))



    ! complex_length in bytes whereas STORAGE_SIZE is in bits
    complex_length = STORAGE_SIZE(coulomb_vertex_I(1)) / 8
    disp_increment = INT(complex_length, KIND=MPI_OFFSET_KIND) * INT(npw, KIND=MPI_OFFSET_KIND)

    call MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(file_in), &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL,unitcv,ierr)


    ! Start with -disp_increment, so that when adding disp_increment, we get 0 in the first iteration
    disp = -disp_increment
    do Ig=1, nI
      disp = disp + disp_increment

      if( ipcol_cd /= INDXG2P(Ig,block_col,0,first_col,npcol_cd) ) cycle
      Il = INDXG2L(Ig,block_col,0,first_col,npcol_cd)

      call MPI_FILE_READ_AT(unitcv, disp, coulomb_vertex_I, &
                            npw, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE,ierr)

      Aread(1:npw,Il)       = coulomb_vertex_I(:)%re
      Aread(npw+1:2*npw,Il) = coulomb_vertex_I(:)%im

      !DEBUG
      !if( Ig == 1 ) then
      !  write(stdout,*) 'Proc', rank,'Testing integral (11|11) (Ha):', DOT_PRODUCT(Aread(:,Il), Aread(:,Il))
      !endif
    enddo

    call MPI_FILE_CLOSE(unitcv, ierr)


    mA = NUMROC(nG, block_row, iprow_sd, first_row, nprow_sd)
    nA = NUMROC(nI, block_col, ipcol_sd, first_col, npcol_sd)
    allocate(A(mA,nA))
    call DESCINIT(descA, nG, nI, block_row, block_col, first_row, first_col, cntxt_sd, mA, info)
    if( info /= 0 ) stop "DESCINIT failure"

    !
    ! Change distribution here
    if( rank == 0 ) write(stdout,'(1x,a,i4,a,i4,a,i4,a,i4,a)') &
                       'Change distribution (', &
                       nprow_cd, ' x ', npcol_cd, ')   to   (', &
                       nprow_sd, ' x ', npcol_sd, ')'
    call PDGEMR2D(nG, nI, Aread, 1, 1, descAread, A, 1, 1, descA, cntxt_sd)

    deallocate(Aread)

    !call PDDOT( nG, rtmp, A, 1, 1, descA, 1, A, 1, 1, descA, 1)
    !DEBUG
    !write(stdout,*) 'Proc', rank, 'Testing integral (11|11) (Ha):', rtmp

    mAt = NUMROC(nI, block_row, iprow_sd, first_row, nprow_sd)
    nAt = NUMROC(nG, block_col, ipcol_sd, first_col, npcol_sd)
    allocate(At(mAt,nAt))
    call DESCINIT(descAt, nI, nG, block_row, block_col, first_row, first_col, cntxt_sd, mAt, info)

    call PDTRAN( nI, nG, 1.0d0, A, 1, 1, descA, 0.0d0, At, 1, 1, descAt )
    deallocate(A)

    !DEBUG
    !call PDDOT( nG, rtmp, At, 1, 1, descAt, nI, At, 1, 1, descAt, nI)
    !write(stdout,*) 'Proc', rank, 'Testing integral (11|11) (Ha):', rtmp

  else
    mAt = NUMROC(nI, block_row, iprow_sd, first_row, nprow_sd)
    nAt = NUMROC(nG, block_col, ipcol_sd, first_col, npcol_sd)
    allocate(At(mAt,nAt))
    call DESCINIT(descAt, nI, nG, block_row, block_col, first_row, first_col, cntxt_sd, mAt, info)

    open(newunit=unitcv, file=TRIM(file_in), form='unformatted', access='stream', status='old', action='read')
    do Ig=1, nI
      read(unitcv) coulomb_vertex_I(:)
        
      At(Ig,1:npw)       = coulomb_vertex_I(:)%re
      At(Ig,npw+1:2*npw) = coulomb_vertex_I(:)%im
    enddo
   
    close(unitcv)
  endif
  call flush(stdout)

end subroutine get_matrix_A


subroutine step1(kp, q, A, descA, Y, descY)
  implicit none

  integer,intent(in) :: kp, q
  real(dp),intent(in)  :: A(:,:)
  integer,intent(in) :: descA(NDEL)
  real(dp),allocatable,intent(out) :: Y(:,:)
  integer,intent(out) :: descY(NDEL)
  !=====
  logical :: read_random
  real(dp) :: start, finish
  integer :: nI, nG
  integer :: iq, info
  integer :: mY, nY
  integer :: mO, nO, descO(NDEL)
  real(dp),allocatable :: Omega(:,:)
  integer :: unitr
  integer :: ikp, ikpl, iGg, iGl
  real(dp) :: rtmp
  !=====

  !
  ! Step 1: Create Y
  !
  if( rank == 0 ) write(*,*) ' **** Step 1 **** '
  call flush(stdout)
  call cpu_time(start)

  nI = descA(M_)
  nG = descA(N_)

  call cpu_time(start)
  ! Omega is nG x kp
  mO = NUMROC(nG, block_row, iprow_sd, first_row, nprow_sd)
  nO = NUMROC(kp , block_col, ipcol_sd, first_col, npcol_sd)
  allocate(Omega(mO,nO))
  call DESCINIT(descO, nG, kp, block_row, block_col, first_row, first_col, cntxt_sd, mO, info)

  ! Y is nI x kp
  mY = NUMROC(nI, block_row, iprow_sd, first_row, nprow_sd)
  nY = NUMROC(kp , block_col, ipcol_sd, first_col, npcol_sd)
  allocate(Y(mY,nY))
  call DESCINIT(descY, nI, kp, block_row, block_col, first_row, first_col, cntxt_sd, mY, info)


  ! Random Omega centered in zero
  inquire(file='random',exist=read_random)
  if( read_random ) then
    if( rank == 0 ) write(stdout,*) 'Read random noise from file', nG, ' x ', kp
    if( rank == 0 ) write(stdout,*) mO, nO
    if( rank == 0 ) write(stdout,*) iprow_sd, nprow_sd, ipcol_sd, npcol_sd
    open(newunit=unitr, file='random', form='unformatted', access='stream', status='old', action='read')
    do ikp=1,kp
      do iGg=1,nG
        read(unitr) rtmp
        if( iprow_sd /= INDXG2P(iGg,block_row,0,first_row,nprow_sd) ) cycle
        if( ipcol_sd /= INDXG2P(ikp,block_col,0,first_col,npcol_sd) ) cycle
        iGl  = INDXG2L(iGg,block_row,0,first_row,nprow_sd)
        ikpl = INDXG2L(ikp,block_col,0,first_col,npcol_sd)
        Omega(iGl,ikpl) = rtmp
      enddo
    enddo
    close(unitr)
    if( rank == 0 ) write(stdout,*) 'Reading done!'
  else
    call random_number(Omega)
    Omega(:,:) = Omega(:,:) - 0.50d0
  endif

  ! Y = A * Omega
  call PDGEMM('N','N', nI, kp, nG, 1.0d0, A, 1, 1, descA, Omega, 1, 1, descO, 0.0d0, Y, 1, 1,descY)
  do iq=1, q
    ! Omega = A**T * Y
    call PDGEMM('T', 'N', nG, kp, nI, 1.0d0, A, 1, 1, descA, Y, 1, 1, descY, 0.0d0, Omega, 1, 1,descO)
    ! Y = A * Omega
    call PDGEMM('N','N', nI, kp, nG, 1.0d0, A, 1, 1, descA, Omega, 1, 1, descO, 0.0d0, Y, 1, 1,descY)
  enddo
  deallocate(Omega)
  call cpu_time(finish)

  if( rank == 0 ) write(stdout,*) 'Step 1: (P)DGEMMs from Omega timing:', finish - start, 'seconds'
  call flush(stdout)

end subroutine step1


subroutine step2(Y, descY, tau)
  implicit none

  real(dp),intent(inout)  :: Y(:,:)
  integer,intent(in)      :: descY(NDEL)
  real(dp),allocatable,intent(out)  :: tau(:)
  !=====
  real(dp) :: start, finish
  integer :: nI, kp
  real(dp),allocatable :: work(:)
  integer :: lwork, info
  !=====

  nI = descY(M_)
  kp = descY(N_)

  call cpu_time(start)
  !
  ! Step 2: Q R decomposition of Y
  !
  if( rank == 0 ) write(*,*) ' **** Step 2 **** '
  if( rank == 0 ) write(*,*) ' Y size:', nI, kp
  call flush(stdout)
  allocate(tau(kp))

  allocate(work(1))
  lwork = -1
  if( nproc == 1 .AND. .FALSE.) then
    call DGEQRF(nI, kp, Y, nI, tau, work, lwork, info)
  else
    call PDGEQRF(nI, kp, Y, 1, 1, descY, tau, work, lwork, info)
  endif
  lwork = INT(work(1))
  deallocate(work)
  if( rank == 0 ) write(*,*) 'allocate workspace:',lwork
  allocate(work(lwork))
  if( rank == 0 ) write(*,*) 'allocate workspace:',ALLOCATED(work)

  if( nproc == 1  .AND. .FALSE.) then
    write(*,*) " DGEQRF before nI, kp", nI, kp
    call DGEQRF(nI, kp, Y, nI, tau, work, lwork, info)
  else
    write(*,*) "PDGEQRF before nI, kp", nI, kp
    call PDGEQRF(nI, kp, Y, 1, 1, descY, tau, work, lwork, info)
  endif
  deallocate(work)

  call cpu_time(finish)
  if( rank == 0 ) write(*,*) 'Step 2: QR timing:', finish - start, 'seconds'
  call flush(stdout)

end subroutine step2

subroutine step3(Y, descY, tau, A, descA, B, descB)
  implicit none

  real(dp),intent(inout)  :: Y(:,:)
  integer,intent(in)      :: descY(NDEL)
  real(dp),allocatable,intent(in)  :: tau(:)
  real(dp),allocatable,intent(inout) :: A(:,:)
  integer,intent(in) :: descA(NDEL)
  real(dp),allocatable :: B(:,:)
  integer,intent(inout) :: descB(NDEL)
  !=====
  real(dp) :: start, finish
  integer :: nI, nG, kp
  real(dp),allocatable :: work(:)
  integer :: lwork, info
  integer :: mB, nB
  !=====
  !

  call cpu_time(start)
  nI = descA(M_)
  nG = descA(N_)
  kp  = descY(N_)


  ! Step 3: Create B
  !
  ! B = Q**T * A
  if( rank == 0 ) write(*,*) ' **** Step 3 **** '
  call flush(stdout)

  ! B is kp x nG
  mB = NUMROC(kp, block_row, iprow_sd, first_row, nprow_sd)
  nB = NUMROC(nG, block_col, ipcol_sd, first_col, npcol_sd)
  allocate(B(mB,nB))
  call DESCINIT(descB, kp, nG, block_row, block_col, first_row, first_col, cntxt_sd, mB, info)

  ! DORMQR applies Q**T on a matrix A
  allocate(work(1))
  lwork = -1
  if( nproc == 1 .AND. .FALSE.) then
    call DORMQR( "L", "T", nI, nG, kp, Y, nI, tau, A, nI, work, lwork, info)
  else
    call PDORMQR( "L", "T", nI, nG, kp, Y, 1, 1, descY, tau, A, 1, 1, descA, work, lwork, info)
  endif
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  if( nproc == 1 .AND. .FALSE.) then
    write(*,*) 'DORMQR'
    call DORMQR( "L", "T", nI, nG, kp, Y, nI, tau, A, nI, work, lwork, info)
  else
    write(*,*) 'PDORMQR'
    call PDORMQR( "L", "T", nI, nG, kp, Y, 1, 1, descY, tau, A, 1, 1, descA, work, lwork, info)
  endif
  deallocate(work)
  if( nproc == 1 .AND. .FALSE.) then
    write(*,*) 'copy'
    B(1:kp,:) = A(1:kp,:)
  else
    !call PDLACPY( " ", kp, nG, A, 1, 1, descA, B, 1, 1, descB)
    write(*,*) 'PDGEMR2D'
    call PDGEMR2D( kp, nG, A, 1, 1, descA, B, 1, 1, descB, cntxt_sd)
  endif

  deallocate(A)

  call cpu_time(finish)
  if( rank == 0 ) write(*,*) 'Step 3: B = Q**T * A product time:', finish - start, 'seconds'
  call flush(stdout)

end subroutine step3


subroutine step4(Y, descY, B, descB, C, descC)
  implicit none

  real(dp),allocatable,intent(inout)  :: Y(:,:)
  integer,intent(in)      :: descY(NDEL)
  real(dp),allocatable,intent(inout)  :: B(:,:)
  integer,intent(in)      :: descB(NDEL)
  real(dp),allocatable,intent(inout)  :: C(:,:)
  integer,intent(out)     :: descC(NDEL)
  !=====
  real(dp) :: start, finish
  integer :: nI, nG, kp, i
  real(dp),allocatable :: sigma(:)
  real(dp),allocatable :: work(:)
  integer :: lwork, info
  integer :: mU, nU
  real(dp),allocatable :: VT(:,:)
  integer :: descVT(NDEL)
  integer :: mC, nC
  real(dp),allocatable :: U(:,:)
  integer :: descU(NDEL)
  !=====

  call cpu_time(start)

  nI = descY(M_)
  nG = descB(N_)
  kp  = descY(N_)
  if( descY(N_) /= descB(M_) ) then
    stop 'Mismatch in matrix dims for Y and B'
  endif

  !
  ! Step 4: SVD of B
  !
  if( rank == 0 ) write(*,*) ' **** Step 4 **** '
  if( rank == 0 ) write(*,*) ' Dimensions nI, k+p, nG:', nI, kp, nG
  call flush(stdout)

  allocate(sigma(kp))

  ! U is kp x kp
  mU = NUMROC(kp , block_row, iprow_sd, first_row, nprow_sd)
  nU = NUMROC(kp , block_col, ipcol_sd, first_col, npcol_sd)
  allocate(U(mU,nU))
  call DESCINIT(descU, kp, kp, block_row, block_col, first_row, first_col, cntxt_sd, mU, info)
  allocate(VT(1,1))
  call DESCINIT(descVT, 1, 1, block_row, block_col, first_row, first_col, cntxt_sd, 1, info)

  allocate(work(1))
  lwork = -1
  if( nproc == 1 .AND. .FALSE.) then
    call DGESVD("S", "N", kp, nG, B, kp, sigma, U, kp, VT, 1, work, lwork, info)
  else
    call PDGESVD("V", "N", kp, nG, B, 1, 1, descB, sigma, U, 1, 1, descU, VT, 1, 1, descVT, work, lwork, info)
  endif
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  if( nproc == 1 .AND. .FALSE.) then
    write(*,*) 'DGESVD'
    call DGESVD("S", "N", kp, nG, B, kp, sigma, U, kp, VT, 1, work, lwork, info)
  else
    if( rank == 0 ) write(*,*) 'PDGESVD'
    call PDGESVD("V", "N", kp, nG, B, 1, 1, descB, sigma, U, 1, 1, descU, VT, 1, 1, descVT, work, lwork, info)
  endif
  deallocate(work)

  !!DEBUG
  !call MPI_BARRIER(MPI_COMM_WORLD,info)
  !write(*,*) "DEBUG"
  !do i=1,kp
  !  write(200+rank,*) i, sigma(i)
  !enddo
  !call flush(200+rank)
  !call MPI_BARRIER(MPI_COMM_WORLD,info)

  if( rank == 0 ) write(*,*) 'PDSCAL'
  do i=1,kp
    call PDSCAL(kp, sigma(i), U, 1, i, descU,1)
    if( rank == 0 ) write(200,*) sigma(i)
  enddo
  call flush(200)

  if( rank == 0 ) write(*,*) 'Singular values:', sigma(1), sigma(kp)
  call cpu_time(finish)
  deallocate(sigma)
  deallocate(B)

  ! Store U in a block of C
  !allocate(C(nI,kp))
  ! C is nI x kp
  mC = NUMROC(nI, block_row, iprow_sd, first_row, nprow_sd)
  nC = NUMROC(kp , block_col, ipcol_sd, first_col, npcol_sd)
  allocate(C(mC,nC))
  call DESCINIT(descC, nI, kp, block_row, block_col, first_row, first_col, cntxt_sd, mC, info)

  C(:,:) = 0.0d0
  if( nproc == 1 .AND. .FALSE.) then
    write(*,*) 'copy'
    C(1:kp,:) = U(1:kp,:)
  else
    write(*,*) 'PDGEMR2D'
    !call PDLACPY( " ", kp, kp, U, 1, 1, descU, C, 1, 1, descC)
    call PDGEMR2D( kp, kp, U, 1, 1, descU, C, 1, 1, descC, cntxt_sd)
  endif
  deallocate(U)
  if( rank == 0 ) write(*,*) 'Step 4: B SVD time:', finish - start, 'seconds'
  call flush(stdout)

end subroutine step4


subroutine step5(Y, descY, tau, C, descC)
  implicit none

  real(dp),allocatable,intent(inout)  :: Y(:,:)
  integer,intent(in)      :: descY(NDEL)
  real(dp),allocatable,intent(inout) :: tau(:)
  real(dp),allocatable,intent(inout)  :: C(:,:)
  integer,intent(out)     :: descC(NDEL)
  !=====
  real(dp) :: start, finish
  integer :: nI, kp, i
  real(dp),allocatable :: sigma(:)
  real(dp),allocatable :: work(:)
  integer :: lwork, info
  !=====

  nI = descY(M_)
  kp  = descY(N_)

  ! Step 5: C = Q * U
  if( rank == 0 ) write(*,*) ' **** Step 5 **** '
  call flush(stdout)
  call cpu_time(start)

  ! DORMQR applies Q on a matrix C
  allocate(work(1))
  lwork = -1
  if( nproc == 1 .AND. .FALSE.) then
    call DORMQR( "L", "N", nI, kp, kp, Y, nI, tau, C, nI, work, lwork, info)
  else
    call PDORMQR( "L", "N", nI, kp, kp, Y, 1, 1, descY, tau, C, 1, 1, descC, work, lwork, info)
  endif
  lwork = INT(work(1))
  deallocate(work)
  if( rank == 0 ) write(*,*) 'lwork:', lwork
  allocate(work(lwork))

  if( nproc == 1 .AND. .FALSE.) then
    call DORMQR( "L", "N", nI, kp, kp, Y, nI, tau, C, nI, work, lwork, info)
  else
    call PDORMQR( "L", "N", nI, kp, kp, Y, 1, 1, descY, tau, C, 1, 1, descC, work, lwork, info)
  endif
  deallocate(work)
  deallocate(tau)
  deallocate(Y)

  call cpu_time(finish)
  if( rank == 0 ) write(*,*) 'Step 5: C = Q * U product time:', finish - start, 'seconds'
  call flush(stdout)

end subroutine step5


subroutine full_svd(A, descA, C, descC)
  implicit none

  real(dp),allocatable,intent(inout)  :: A(:,:)
  integer,intent(in)      :: descA(NDEL)
  real(dp),allocatable,intent(inout)  :: C(:,:)
  integer,intent(out)     :: descC(NDEL)
  !=====
  real(dp) :: start, finish
  integer :: nI, nG, i
  real(dp),allocatable :: sigma(:)
  real(dp),allocatable :: work(:)
  integer :: lwork, info
  real(dp),allocatable :: VT(:,:)
  integer :: descVT(NDEL)
  integer :: mC, nC
  !=====

  call cpu_time(start)

  nI = descA(M_)
  nG = descA(N_)

  !
  ! Direct SVD of C
  !
  if( rank == 0 ) write(*,*) ' **** Direct SVD **** '

  allocate(sigma(nG))

  ! C is nI x nG
  mC = NUMROC(nI, block_row, iprow_sd, first_row, nprow_sd)
  nC = NUMROC(nG, block_col, ipcol_sd, first_col, npcol_sd)
  allocate(C(mC,nC))
  allocate(VT(1,1))
  call DESCINIT(descC, nI, nG, block_row, block_col, first_row, first_col, cntxt_sd, mC, info)

  allocate(work(1))
  lwork = -1
  call PDGESVD("V", "N", nI, nG, A, 1, 1, descA, sigma, C, 1, 1, descC, VT, 1, 1, descVT, work, lwork, info)
  lwork = INT(work(1))
  deallocate(work)
  allocate(work(lwork))

  call PDGESVD("V", "N", nI, nG, A, 1, 1, descA, sigma, C, 1, 1, descC, VT, 1, 1, descVT, work, lwork, info)
  deallocate(work)

  !!DEBUG
  !call MPI_BARRIER(MPI_COMM_WORLD,info)
  !write(*,*) "DEBUG"
  !do i=1,k
  !  write(200+rank,*) i, sigma(i)
  !enddo
  !call flush(200+rank)
  !call MPI_BARRIER(MPI_COMM_WORLD,info)

  do i=1,nG
    !C(:,i) = C(:,i) * sigma(i)
    call PDSCAL(nG, sigma(i), C, 1, i, descC,1)
    if( rank == 0 ) write(200,*) sigma(i)
  enddo
  if( rank == 0 ) write(*,*) 'Singular values:', sigma(1), sigma(nG)
  call cpu_time(finish)
  deallocate(sigma)
  deallocate(A)

  if( rank == 0 ) write(*,*) 'A SVD time:', finish - start, 'seconds'

end subroutine full_svd

subroutine dump_matrix_C(k, file_out, C, descC)
  implicit none

  integer,intent(in) :: k
  character(len=*),intent(in) :: file_out
  real(dp),allocatable,intent(inout)  :: C(:,:)
  integer,intent(out)     :: descC(NDEL)
  !=====
  integer :: nI, kp
  real(dp),allocatable :: Ct(:,:)
  integer :: descCt(NDEL)
  integer :: mCt, nCt
  real(dp),allocatable :: Ctdump(:,:)
  integer :: descCtdump(NDEL)
  integer :: mCtdump, nCtdump
  integer(kind=MPI_OFFSET_KIND) :: disp, disp_increment
  integer :: unitcv, ierr, info
  complex(dp),allocatable :: coulomb_vertex_I(:)
  integer :: complex_length, kc, Ig, Il
  !=====

  nI = descC(M_)
  kp = descC(N_)
  kc = k / 2

  if( rank == 0 ) write(stdout,*) 'C global matrix:', nI, ' x ', kp
  if( rank == 0 ) write(stdout,*) 'C local matrix:', SIZE(C,DIM=1), ' x ', SIZE(C,DIM=2)
  call flush(stdout)


  mCt = NUMROC(kp, block_row, iprow_sd, first_row, nprow_sd)
  nCt = NUMROC(nI, block_col, ipcol_sd, first_col, npcol_sd)
  allocate(Ct(mCt,nCt))
  call DESCINIT(descCt, kp, nI, block_row, block_col, first_row, first_col, cntxt_sd, mCt, info)
  if( info /= 0 ) stop "DESCINIT descCt failure"

  ! Ct = C**T
  call PDTRAN( kp, nI, 1.0d0, C, 1, 1, descC, 0.0d0, Ct, 1, 1, descCt )
  deallocate(C)

  mCtdump = NUMROC(k , block_row, iprow_cd, first_row, nprow_cd)
  nCtdump = NUMROC(nI, block_col, ipcol_cd, first_col, npcol_cd)
  allocate(Ctdump(mCtdump,nCtdump))
  call DESCINIT(descCtdump, k, nI, block_row, block_col, first_row, first_col, cntxt_cd, mCtdump, info)
  if( info /= 0 ) stop "DESCINIT descCtdump failure"

  !
  ! Change distribution here
  if( rank == 0 ) write(stdout,'(1x,a,i4,a,i4,a,i4,a,i4,a)') &
                     'Change distribution (', &
                     nprow_sd, ' x ', npcol_sd, ')   to   (', &
                     nprow_cd, ' x ', npcol_cd, ')'
  if( nproc == 1 ) then
    Ctdump(1:k,:) = Ct(1:k,:)
  else
    call PDGEMR2D( k, nI, Ct, 1, 1, descCt, Ctdump, 1, 1, descCtdump, cntxt_sd)
  endif
  deallocate(Ct)


  if( rank == 0 ) write(stdout,'(/,1x,a,a,a)') 'Writing file ', TRIM(file_out), ' with MPI-IO'
  if( rank == 0 ) write(stdout,'(5x,a,i4,a,i4)') 'using a processor grid:', nprow_cd, ' x ', npcol_cd
  if( rank == 0 ) write(stdout,'(5x,a,i6,a,i8)') 'CoulombVertex dimensions:', kc, ' x ', nI

  allocate(coulomb_vertex_I(kc))

  ! complex_length in bytes whereas STORAGE_SIZE is in bits
  complex_length = STORAGE_SIZE(coulomb_vertex_I(1)) / 8
  disp_increment = INT(complex_length, KIND=MPI_OFFSET_KIND) * INT(kc, KIND=MPI_OFFSET_KIND)

  ! Erase file first
  call MPI_FILE_DELETE(TRIM(file_out), info, ierr)
  if( rank == 0 ) write(stdout,*) "File size (bytes):",INT(nI,KIND=MPI_OFFSET_KIND) * disp_increment
  call MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(file_out), &
                     MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                     MPI_INFO_NULL,unitcv,ierr)
  if( ierr /= 0 ) stop 'error opening file'

  ! Start with -disp_increment, so that when adding disp_increment, we get 0 in the first iteration
  disp = -disp_increment
  do Ig=1, nI
    disp = disp + disp_increment

    if( ipcol_cd /= INDXG2P(Ig,block_col,0,first_col,npcol_cd) ) cycle
    Il = INDXG2L(Ig,block_col,0,first_col,npcol_cd)

    coulomb_vertex_I(:) = CMPLX(Ctdump(1:kc,Il), Ctdump(kc+1:2*kc,Il))

    call MPI_FILE_WRITE_AT(unitcv, disp, coulomb_vertex_I, &
                          kc, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE,ierr)


  enddo
  call MPI_FILE_CLOSE(unitcv, ierr)
  call flush(stdout)

end subroutine dump_matrix_C

end module low_level


