module input_param
  use low_level

  integer, protected :: p, k, q
  integer, protected :: nI, nG, npw
  character(len=128), protected :: file_in, file_out, method
  namelist /input/ nI, nG, npw, k, q, p, file_in, file_out, method
  
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

  if( COMMAND_ARGUMENT_COUNT() == 1 ) then
    call GET_COMMAND_ARGUMENT(1, VALUE=input_file_name)
  else
    stop "no input file"
  endif

  if( rank == 0 ) write(stdout,*) 'Reading input file: ', input_file_name
  open(newunit=ifile, file=TRIM(input_file_name), status='old', action='read')
  read(ifile,input)
  close(ifile)

  if( npw == 0 .AND. nG == 0 ) stop 'Set npw or nG'
  if( npw /= 0 .AND. nG /= 0 )  stop 'Set npw or nG, not both'
  if( nG  /= 0 ) npw = nG / 2
  if( npw /= 0 ) nG = npw * 2

  ! Assume nI >> nG
  ! Enforce it
  if( nI <= nG ) stop "nI <= nG"

  if( rank == 0 ) write(*,*) 'k=',k
  if( rank == 0 ) write(*,*) 'p=',p
  if( rank == 0 ) write(*,*) 'q=',q

  kp = k + p

end subroutine read_input_file

end module input_param
