subroutine print_hello(mpirank)
use iso_fortran_env
!! use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
!!                                           stdout=>output_unit, &
!!                                           stderr=>error_unit

implicit none

integer mpirank, stdout
character(len=*), parameter :: nullfile = '/dev/null'

if (mpirank > 0) then
	open(output_unit,file=nullfile)
endif

write(*,'(A,I5)') 'Hello from rank ', mpirank


end