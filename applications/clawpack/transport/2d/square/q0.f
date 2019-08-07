c     # ---------------------------------------------------------------      
      double precision function q0_physical(xp,yp,zp)
      implicit none

      double precision xp, yp, zp

      double precision x0, y0, z0, q0, r0, r
      double precision xpp, ypp
      double precision Hsmooth

c     # Sphere centered at (0.5,0.5,0) on swirl
      x0 = 0.5
      y0 = 0.5
      z0 = 0
      r0 = 0.2   !! radius of sphere

c     # modulo(A,P) = A - floor(A/P)*P
c      write(6,*) 'mod(6.1,1.d0)  ', modulo(6.1d0,1.d0)
c      write(6,*) 'mod(0.5,1.d0)  ', modulo(0.5d0,1.d0)
c      write(6,*) 'mod(-0.5,1.d0) ', modulo(-0.5d0,1.d0)
c      write(6,*) 'mod(-1.3,1.d0) ', modulo(-1.3d0,1.d0)
c      write(6,*) ' '
c      write(6,*) 'mod(6.1,1.d0)  ', 6.1d0 - floor(6.1d0)
c      write(6,*) 'mod(0.5,1.d0)  ', 0.5d0 - floor(0.5d0)
c      write(6,*) 'mod(-0.5,1.d0) ', -0.5d0 - floor(-0.5d0)
c      write(6,*) 'mod(-1.3,1.d0) ', -1.3d0 - floor(-1.3d0)
c      stop

      xpp = modulo(xp,1.d0)
      ypp = modulo(yp,1.d0)


      r = sqrt((xpp - x0)**2 + (ypp-y0)**2)
      q0 = Hsmooth(r + r0) - Hsmooth(r - r0)

      q0_physical = q0

      end

c     # ---------------------------------------------------------------      
      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end



