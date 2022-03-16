c
c
c
c     =================================================
      function fdisc(blockno,xc,yc)
c     =================================================
      implicit double precision (a-h,o-z)
      common/cdisc/ x0,y0,alf,beta,r0,idisc


      integer blockno
      !integer*8 cont  get_context
      double precision f1,f2

c
c     # for computing cell averages for initial data that has a
c     # discontinuity along some curve.  fdisc should be negative to the
c     # left of the curve and positive to the right
c     # idisc specifies the nature of the discontinuity for two
c     # particular cases (a straight line and circle) but this routine
c     # can be modified for any other curve.
c

      x = xc
      y = yc

      go to (10,20,30,40) idisc
c
   10 continue
c     # straight line through (x0,y0) with normal (alf,beta) pointing
c     # into right state
c
      fdisc = (x-x0)*alf + (y-y0)*beta
      return
c
   20 continue
c     # circle of radius r0:
      fdisc = (x-x0)**2 + (y-y0)**2 - r0**2
      return

   30 continue
      f1 = (x-x0)**2 + (y-y0+0.5d0)**2 - r0**2
      f2 = (x-x0)**2 + (y-y0-0.5d0)**2 - r0**2
      if (f1 .lt. 0) then
         fdisc = f1
      elseif (f2 .lt. 0) then
         fdisc = f2
      else
         fdisc = min(f1,f2)
      endif
      return

   40 continue
      f1 = (x-x0)**2 + (y-y0+0.5d0)**2 - (r0/sqrt(2.d0))**2
      f2 = (x-x0)**2 + (y-y0-0.5d0)**2 - (r0/sqrt(2.d0))**2
      if (f1 .lt. 0) then
         fdisc = f1
      elseif (f2 .lt. 0) then
         fdisc = f2
      else
         fdisc = min(f1,f2)
      endif


      return
      end
