      subroutine setprob
      implicit double precision (a-h,o-z)
      character(len=25) fname

      common /compsi/ pi
      common /comvt/ tperiod,pi2

c
c     # compute pi, used in psi.f
      pi = 4.d0 * datan(1.d0)
c
c     # save 2*pi and tperiod in common block for use in b4step2:
c
      pi2 = 2.d0*pi
c     # Add 2016/4/18 tperiod
      tperiod = 4.d0
c
      iunit = 7
c      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
c      call opendatafile(iunit, fname)

c      read(7,*) tperiod

      return
      end
