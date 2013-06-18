c=========================================================================
      subroutine setgauges
c=========================================================================

      implicit double precision (a-h,o-z)
      character*20 fname
      logical foundFile

      include "gauges.i"
      include "call.i"


      fname  = 'setgauges.data'
      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname 
        stop
      endif

      open(unit=7,file=fname,status='old',form='formatted')

      read(7,*) mgauges
      if (mgauges.gt.maxgauges) then
            write(*,*) 'ERROR in setgauges'
            write(*,*) 'mgauges = ',mgauges,'   maxgauges = ',maxgauges
            write(*,*) 'increase maxgauges in gauges.i'
            stop
            endif

      do i=1,mgauges
          read(7,*) igauge(i),xgauge(i),ygauge(i)
          mbestsrc(i) = 0   ! initialize for starters
          enddo
      close(7)

c     # open file for output of gauge data
c     # all data is output in one binary file with format
c     # gauge number, time, depth

      open(unit=OUTGAUGEUNIT,file='fort.gauge',status='unknown',
     .           form='formatted')

      return
      end
