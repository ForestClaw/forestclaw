subroutine print_hello(mpirank, verb)
use iso_fortran_env
implicit none

integer mpirank, i, verb
character(len=*), parameter :: nullfile = '/dev/null'

if (mpirank > 0) then
	open(output_unit,file=nullfile)
endif

write(*,'(A,I5)') "Hello from rank ", mpirank

do i = 1,10
	write(*,'(I5)') i
enddo

if (mpirank > 0) then
	close(output_unit)
endif

end