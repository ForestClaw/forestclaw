      subroutine mapc2m(xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp

      double precision ll(2), ur(2)
      logical isperiodic


c     # This constructs mapping in [-1,-1]x[1,1].  To get
c     # something in a box of a different size and location
c     # translate below.

      isperiodic = .false.
      call mapc2m_inclusions(xc,yc,xp,yp,isperiodic)

c     # Map to from [-1,1]x[-1,1] to [ll]x[ur]

c     # The transformation below doesn't do anything, but demonstrates
c     # how to translate and scale the output from the inclusions
c     # maping.
      ll(1) = -1
      ll(2) = -1
      ur(1) = 1
      ur(2) = 1
      call transform(xp,yp,ll,ur)

      zp = 0


      end

      subroutine transform(xp,yp,ll,ur)
      implicit none

      double precision xp,yp,ll(2),ur(2)

      xp = ll(1) + (ur(1) - ll(1))*(xp + 1.d0)/2.d0
      yp = ll(2) + (ur(2) - ll(2))*(yp + 1.d0)/2.d0

      end

      double precision function exact_area()
      implicit none

      exact_area = 4.d0
      end
