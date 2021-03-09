      subroutine wavetank_setprob()
      implicit none

      DOUBLE PRECISION pi, pi2
      COMMON /compi/ pi, pi2

      DOUBLE PRECISION  grav, dry_tolerance, sea_level
      COMMON /common_swe/ grav, dry_tolerance, sea_level

      DOUBLE PRECISION breaking, alpha, g_sgn
      COMMON /common_sgn/ breaking, alpha, g_sgn

c     # These should be read in as options.

      pi = 4.d0*atan(1.d0)
      
      open(10,file='setprob.data')
      read(10,*) grav
      read(10,*) dry_tolerance
      read(10,*) sea_level

c     # SGN parameters      
      read(10,*) breaking
      read(10,*) alpha
      close(10)

      g_sgn = grav

      return
      end
