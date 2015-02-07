      subroutine setprob
      implicit double precision (a-h,o-z)
      common /comaux/ rhol,cl,rhor,cr

c     # this should be set from user options.
      rhol = 1.d0
      cl = 1.d0
      rhor = 4.d0
      cr = 0.5d0

      return
      end
