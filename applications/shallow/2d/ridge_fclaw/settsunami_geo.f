c=========================================================================
      subroutine settsunami
c=========================================================================

      implicit double precision (a-h,o-z)
      character*20 fname
      logical foundFile

      include "call.i"
      include "geo.i"

      write(parmunit,*) ' '
      write(parmunit,*) '--------------------------------------------'
      write(parmunit,*) 'SETTSUNAMI:'
      write(parmunit,*) '-----------'

c       # read user parameters from settsunami.data

      fname  = 'settsunami.data'
      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname
        stop
      endif

      open(unit=7,file=fname,status='old',form='formatted')

      read(7,*) sealevel
      read(7,*) drytolerance
      read(7,*) wavetolerance
      read(7,*) depthdeep
      read(7,*) maxleveldeep
      read(7,*) ifriction
      read(7,*) coeffmanning
      read(7,*) frictiondepth
      close(7)

      write(parmunit,*) '   drytolerance:',drytolerance
      write(parmunit,*) '   wavetolerance:',wavetolerance
      write(parmunit,*) '   maxleveldeep:', maxleveldeep
      write(parmunit,*) '   depthdeep:', depthdeep
      write(parmunit,*) '   ifriction:', ifriction
      write(parmunit,*) '   Manning coefficient:',coeffmanning
      write(parmunit,*) '   frictiondepth:',frictiondepth

      return
      end
