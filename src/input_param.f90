module input_param
  use low_level

  integer, protected :: p, k, q
  integer, protected :: nI, nG, npw, nmo, nmo_file
  logical, protected :: just_check
  character(len=128), protected :: file_in, file_out, method
  namelist /input/ nI, nG, nmo, nmo_file, npw, k, q, p, file_in, file_out, method, just_check
  
  integer, protected :: kp

contains

subroutine read_input_file()
  implicit none
  integer :: ifile
  character(len=128) :: input_file_name

  ! Default
  method = 'PROJ_QR_SVD'
  k = 1000
  p = 0
  q = 1
  npw = 0
  nG  = 0
  nmo = 0
  nI  = 0
  nmo_file = 0
  just_check = .FALSE.

  if( COMMAND_ARGUMENT_COUNT() == 1 ) then
    call GET_COMMAND_ARGUMENT(1, VALUE=input_file_name)
  else
    stop "no input file"
  endif

  if( rank == 0 ) write(stdout, *) 'Reading input file: ', input_file_name
  open(newunit=ifile, file=TRIM(input_file_name), status='old', action='read')
  read(ifile, input)
  close(ifile)

  if( npw == 0 .AND. nG == 0 ) stop 'Set npw or nG'
  if( npw /= 0 .AND. nG /= 0 )  stop 'Set npw or nG, not both'
  if( nG  /= 0 ) npw = nG / 2
  if( npw /= 0 ) nG = npw * 2

  if( nmo == 0 .AND. nI == 0 ) stop 'Set nmo or nI'
  if( nmo /= 0 .AND. nI /= 0 ) stop 'Set nmo or nI, not both'
  if( nmo /= 0 ) nI = nmo**2
  if( nI /= 0 )  nmo = INT( SQRT( REAL(nI) ) )
  if( nmo_file == 0 ) nmo_file = nmo

  ! Assume nI >> nG
  ! Enforce it
  if( nI <= nG ) then
    write(*, *) nmo, nI, nG
    stop "Only for rectangular matrices with more rows than columns: nI <= nG"
  endif

  if( rank == 0 ) write(*, *) 'k=', k
  if( rank == 0 ) write(*, *) 'p=', p
  if( rank == 0 ) write(*, *) 'q=', q

  kp = k + p

end subroutine read_input_file

end module input_param
