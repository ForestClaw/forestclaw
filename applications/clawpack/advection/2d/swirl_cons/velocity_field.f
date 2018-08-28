      subroutine velocity_field(xc,yc,u,v)
      implicit none

      double precision xc,yc,u,v

      double precision s, pi
      double precision xp,yp,zp
      integer example

      common /compi/ pi
      common /comex/ example

      s = sqrt(2.d0)
      if (example .eq. 0) then
c        # Conservative for all solvers (rp=1,2,3,4)              
c         u = cos(2*pi*xc) + 2
c         v = cos(2*pi*yc) + 2
         u = 1
         v = 1
      elseif (example .eq. 1) then
c        # Conservative for all solvers (rp=1,2,3,4)               
         u = s*(cos(pi*xc)**2 + 0.5d0)
         v = s*(sin(pi*yc)**2 + 0.5d0)
      elseif (example .eq. 2) then
         u = s*(cos(pi*xc)**2 - 0.5d0)
         v = s*(sin(pi*yc)**2 - 0.5d0)
      else if (example .ge. 3) then
c         u = 1.d0
        u = s*(cos(pi*xc)**2 + 0.5d0)

c         v = 0.0d0
         v = s*(sin(pi*yc)**2 + 0.5d0)

      else
         write(6,'(A,A)') 'clawpack46_setaux : ',
     &              'No valid example provided'
         stop
      endif



      end
