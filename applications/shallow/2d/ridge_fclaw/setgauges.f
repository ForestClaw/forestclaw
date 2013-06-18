c=========================================================================
      subroutine setgauges
c=========================================================================

      implicit double precision (a-h,o-z)
      character*20 fname
      logical foundFile

c      include "gauges.i"
c      include "call.i"
c
c
c      fname  = 'setgauges.data'
c      inquire(file=fname,exist=foundFile)
c      if (.not. foundFile) then
c        write(*,*) 'You must provide a file ', fname
c        stop
c      endif
c
c      open(unit=7,file=fname,status='old',form='formatted')
c
c      read(7,*) mgauges
c      if (mgauges.gt.maxgauges) then
c            write(*,*) 'ERROR in setgauges'
c            write(*,*) 'mgauges = ',mgauges,'   maxgauges = ',maxgauges
c            write(*,*) 'increase maxgauges in gauges.i'
c            stop
c            endif
c
c      do i=1,mgauges
c          read(7,*) igauge(i),xgauge(i),ygauge(i)
c          mbestsrc(i) = 0   ! initialize for starters
c          enddo
c      close(7)
c
cc     # open file for output of gauge data
cc     # all data is output in one binary file with format
cc     # gauge number, time, depth
c
c      open(unit=OUTGAUGEUNIT,file='fort.gauge',status='unknown',
c     .           form='formatted')
c
      return
      end
