      subroutine velocity_field(xc,yc,u,v)
      implicit none

      double precision xc,yc,u,v

      double precision s, pi
      double precision xp,yp,zp, r, th
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
c         u = 1
c         u = s*(sin(pi*yc)**2 + 0.5d0)
c         u = s*sin(pi*yc)

c         v = 0
c         v = -s*(sin(pi*xc)**2 + 0.5d0)
c         v = s*sin(pi*xc)

c         r = sqrt((x-0.5d0)**2 + (y-0.5d0)**2)

c        # Rigid body rotation
c        psi = r**2
c         u = 2*(yc-0.5)
c         v = -2*(xc-0.5)
 
c        # Filament formation (negative for clockwise rotation)
c         psi = (4.d0/3.d0)*r**3

      else
         write(6,'(A,A)') 'clawpack46_setaux : ',
     &              'No valid example provided'
         stop
      endif



      end
