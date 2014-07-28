      subroutine mapc2m_inclusions(xc,yc,xp,yp,isperiodic)
      implicit none

      double precision xc,yc,xp,yp
      logical isperiodic

      double precision xc1,yc1,xp1,yp1
      double precision xlow,ylow,xhi,yhi, r1
      double precision xc0, yc0
      integer i, nbox
      logical inbox

      double precision boxes_com(6,100)
      integer n_box_com

      common /combox0/ n_box_com
      common /combox1/ boxes_com

      nbox = n_box_com

      xc0 = xc
      yc0 = yc

      if (isperiodic) then
         if (xc .lt. -1.d0) then
            xc0 = xc + 2
         elseif (xc .gt. 1.d0) then
            xc0 = xc - 2
         endif

         if (yc .lt. -1.d0) then
            yc0 = yc + 2
         elseif (yc .gt. 1.d0) then
            yc0 = yc - 2
         endif
      endif

      xp = xc0
      yp = yc0


      inbox = .false.
      do i = 1,nbox
c        # [xlow,ylow]x[xhi,yhi] are the ll and ur corners of
c        # included box in [-1,1]x[-1,1] units.
         xlow = boxes_com(1,i)
         ylow = boxes_com(2,i)
         xhi = boxes_com(3,i)
         yhi = boxes_com(4,i)

c        # Radius of enclosed circle is in 0 < r1 < 1.
         r1 = boxes_com(5,i)

c        # Scale [xlow,ylow]x[xhi,yhi] into [-1,1]x[-1,1]
         xc1 = -1.d0 + 2.d0*(xc0 - xlow)/(xhi - xlow)
         yc1 = -1.d0 + 2.d0*(yc0 - ylow)/(yhi - ylow)

         if (abs(xc1) .le. 1.d0 .and. abs(yc1) .le. 1.d0) then
            call mapc2m_inclusion2(xc1,yc1,xp1,yp1,r1)

c           # Scale physical coordinates back so that mapped
c           # region is in [xlow,ylow]x[xhi,yhi]
            xp = (xp1 + 1.d0)*(xhi - xlow)/2.d0 + xlow
            yp = (yp1 + 1.d0)*(yhi - ylow)/2.d0 + ylow
            inbox = .true.
         endif
         if (inbox) then
            exit
         endif
      enddo
      if (isperiodic) then
         if (xc .lt. -1) then
            xp = xp - 2
         elseif (xc .gt. 1) then
            xp = xp + 2
         endif
         if (yc .lt. -1) then
            yp = yp - 2
         elseif (yc .gt. 1) then
            yp = yp + 2
         endif
      endif

      end

c     # ------------------------------------------------------------------
c     MAPC2M_INCLUSION maps the grid [-1,1]x[-1,1] to a grid enclosing a
c     circle of radius r1 < 1.
c     # ------------------------------------------------------------------

      subroutine mapc2m_inclusion2(xc,yc,xp,yp,r1)
      implicit none

      double precision xc,yc,xp,yp,zp,r1

      double precision d, dinv, scale, s, r, w, ct
      logical use_convex

      d = max(abs(xc),abs(yc))
      r = max(sqrt(xc**2 + yc**2),1e-10)

      if (d .lt. r1) then
c        # --------------------------------------------------------------
c        # Map area inside of the circle
c        # -------------------------------------------------------------------

         scale = d/r

         xp = scale*xc
         yp = scale*yc

c        # Blend this map with the Cartesian map to get
c        # more Cartesian like grid near the center of the
c        # circle.  If we set this to false, we get concentric
c        # circles.
         use_convex = .true.
         if (use_convex) then
            w = (d/r1)**2
            xp = w*xp + (1.d0-w)*xc/sqrt(2.d0)
            yp = w*yp + (1.d0-w)*yc/sqrt(2.d0)
         endif
      else

c        # --------------------------------------------------------------
c        # Map area outside of the circle
c        # -------------------------------------------------------------------
         dinv = 1.d0/d

c        # d/r = cos(theta), where theta is the angle that a point
c        # makes with the x-axis.

         ct = d/r
         s = (1.d0 - r1*ct)/(1.d0 - r1)
         scale = s*(1.d0 - dinv) + dinv

         xp = scale*xc
         yp = scale*yc

      endif

      end
