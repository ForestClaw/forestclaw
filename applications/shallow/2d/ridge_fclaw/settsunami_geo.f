c=========================================================================
      subroutine settsunami
c=========================================================================

      implicit double precision (a-h,o-z)
      character*20 fname
      logical foundFile

c      include "call.i"
c      include "geo.i"
c
c      write(parmunit,*) ' '
c      write(parmunit,*) '--------------------------------------------'
c      write(parmunit,*) 'SETTSUNAMI:'
c      write(parmunit,*) '-----------'
c
cc       # read user parameters from settsunami.data
c
      fname  = 'settsunami.data'
c      inquire(file=fname,exist=foundFile)
c      if (.not. foundFile) then
c        write(*,*) 'You must provide a file ', fname
c        stop
c      endif

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

c      write(parmunit,*) '   drytolerance:',drytolerance
c      write(parmunit,*) '   wavetolerance:',wavetolerance
c      write(parmunit,*) '   maxleveldeep:', maxleveldeep
c      write(parmunit,*) '   depthdeep:', depthdeep
c      write(parmunit,*) '   ifriction:', ifriction
c      write(parmunit,*) '   Manning coefficient:',coeffmanning
c      write(parmunit,*) '   frictiondepth:',frictiondepth

      return
      end
