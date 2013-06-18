      subroutine setprob
      implicit double precision (a-h,o-z)
      common /comsphere/ Rsphere, Omega
      common /comic/ Px,Py,Pz
      common /cprob/ ampl, tol1, tol2
      common /sw/  g

      double precision rot_angle(2)

      common /compi/ pi

      pi = 4.d0*atan(1.d0)

      call settsunami
      call setgauges
c
      open(unit=7,file='setprob.data',status='old',form='formatted')
c
c     # radius of sphere:
      Rsphere = 6.371d6

c     # rotation rate:
      Omega = 7.292d-5 ! 1/sec
c     Omega = 0.0d0    ! 1/sec
c
c     # graviational constant
c     g =  11489.57219e0
      g =  1.d0  ! when using geopotential height

      read(7,*) x1, y1, z1
c     # normalize to have length 1 in case it didn't:
      Px = x1
      Py = y1
      Pz = z1
      Pnorm = dsqrt(Px**2 + Py**2 + Pz**2)
      Px = Px / Pnorm
      Py = Py / Pnorm
      Pz = Pz / Pnorm

c     # amplitude of initial Gaussian:
      read(7,*) ampl
c     # refinement tolerances used in flag2refine:
      read(7,*) tol1, tol2

c     # Set mapping scaling
      scale = Rsphere
      rot_angle(1) = 0
      rot_angle(2) = 0
      call setup_mappedgrid(rot_angle,scale)


      return
      end
