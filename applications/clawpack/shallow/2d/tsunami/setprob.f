      subroutine tsunami_setprob()
      implicit none

      DOUBLE PRECISION pi, pi2
      COMMON /compi/ pi, pi2

      DOUBLE PRECISION a,b,h0
      COMMON /common_initheight/ a,b,h0

      DOUBLE PRECISION  grav, dry_tolerance, sea_level
      COMMON /common_swe/ grav, dry_tolerance, sea_level

      DOUBLE PRECISION breaking, alpha, g_sgn
      COMMON /common_sgn/ breaking, alpha, g_sgn

c     # These should be read in as options.
      
      open(10,file='setprob.data')

      !! User options
      read(10,*) grav
      read(10,*) a
      read(10,*) b
      read(10,*) h0
      read(10,*) dry_tolerance
      read(10,*) sea_level

      !! From SGN options
      read(10,*) breaking
      read(10,*) alpha
      close(10)

      g_sgn = grav

      return
      end
