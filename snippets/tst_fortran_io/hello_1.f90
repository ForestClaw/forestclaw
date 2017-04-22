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
integer fclaw_debug
character(len=*), parameter :: nullfile = '/dev/null'

silent = 0
essential = 1
production = 2
info = 3
debug = 4

fclaw_global_essential = 50
fclaw_global_production = 51
fclaw_global_info = 52
fclaw_debug = 53

if (verb .eq. debug) then
	fclaw_debug = output_unit
	if (mpirank .eq. 0) then
		fclaw_global_info = output_unit
		fclaw_global_production = output_unit
		fclaw_global_essential = output_unit
	else
		!! ranks > 0 don't print at all
		open(fclaw_global_info,file=nullfile)			
		open(fclaw_global_production,file=nullfile)			
		open(fclaw_global_essential,file=nullfile)			
	endif
elseif (verb .eq. info) then
	open(fclaw_debug,file=nullfile)			
	if (mpirank .eq. 0) then
		fclaw_global_info = output_unit
		fclaw_global_production = output_unit
		fclaw_global_essential = output_unit
	else
		open(fclaw_global_info,file=nullfile)			
		open(fclaw_global_production,file=nullfile)			
		open(fclaw_global_essential,file=nullfile)			
	endif
elseif (verb .eq. production) then		
	open(fclaw_debug,      file=nullfile)			
	open(fclaw_global_info,file=nullfile)
	if (mpirank == 0) then
		fclaw_global_production = output_unit
		fclaw_global_essential = output_unit
	else			
		open(fclaw_global_production,file=nullfile)			
		open(fclaw_global_essential,file=nullfile)			
	endif
elseif (verb .eq. essential) then
	open(fclaw_debug,            file=nullfile)	
	open(fclaw_global_info,      file=nullfile)
	open(fclaw_global_production,file=nullfile)
	if (mpirank .eq. 0) then
		fclaw_global_essential = output_unit
	else
		open(fclaw_global_essential,file=nullfile)			
	endif
elseif (verb .eq. silent) then
	open(output_unit,file=nullfile)
endif

!! write(*,'(A,I5)') '(*) Hello from rank ', mpirank
write(fclaw_global_silent,*)     "(silent) Hello World!"
write(fclaw_global_essential,*)  "(essential) Hello World! "
write(fclaw_global_production,*) "(production) Hello World! "
write(fclaw_global_info,*)       "(info) Hello World! "
write(fclaw_debug,*)             "(debug) Hello from rank ", mpirank

close(fclaw_global_essential)
close(fclaw_global_production)
close(fclaw_global_info)
close(fclaw_debug)


close(output_unit)

end