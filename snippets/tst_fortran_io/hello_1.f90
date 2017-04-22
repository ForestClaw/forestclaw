subroutine print_hello_1(mpirank,verb)
use iso_fortran_env
!! use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
!!                                           stdout=>output_unit, &
!!                                           stderr=>error_unit

implicit none

integer mpirank, verb
integer silent, essential, production, info, debug
integer fclaw_global_info, fclaw_global_production
integer fclaw_global_essential, fclaw_global_silent
integer fclaw_debug, nullfile_unit, r

character(len=*), parameter :: nullfile = '/dev/null'
character(len=*), parameter :: fmt_info       = '("[ash3d] ",A)'
character(len=*), parameter :: fmt_production = '("[ash3d] ",A)'
character(len=*), parameter :: fmt_essential  = '("[ash3d] ",A)'
character(len=*), parameter :: fmt_debug      = '("[",I5,"] ",A)'

silent = 0
essential = 1
production = 2
info = 3
debug = 4

nullfile_unit = 50

!!mpirank = mpirank + 30045

!!k = 1
!!r = mpirank
!!do while (.true.)
!!	r = r/10
!!	write(6,*) r
!!	k = k + 1
!!	if (r .eq. 0) then
!!		exit
!!	endif
!!enddo	
!!write(6,*) 'Digits = ', k

if (verb .gt. 0) then
	open(nullfile_unit,file=nullfile)	
else
	open(output_unit,file=nullfile)
endif

fclaw_global_essential = nullfile_unit
fclaw_global_production = nullfile_unit
fclaw_global_info = nullfile_unit
fclaw_debug = nullfile_unit

!! Every processor writes something in debug mode
if (verb .eq. debug) then
	fclaw_debug  = output_unit
endif

!! Only processor 0 writes in other modes
if (mpirank == 0) then
	if (verb .ge. essential) then
		fclaw_global_essential = output_unit
		if (verb .ge. production) then
			fclaw_global_production = output_unit
			if (verb .ge. info) then
				fclaw_global_info = output_unit
			endif
		endif
	endif
endif

write(*,'(A,I5)')                "(stdout) Hello, World! from rank ", mpirank


write(fclaw_global_silent,*)                  "(silent)     Hello World!"
write(fclaw_global_essential,fmt_essential)   "(essential)  Hello World! "
write(fclaw_global_production,fmt_production) "(production) Hello World! "
write(fclaw_global_info,fmt_info)             "(info)       Hello World! "
write(fclaw_debug,fmt_debug)                  mpirank+200,      "Hello, World!"

end