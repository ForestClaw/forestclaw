c -*- FORTRAN -*-
c     Common blocks needed by step2, flux2.  We have these here in a common
C     block, rather than as arguments, so that we can use the same step2.f and
C     flux2.f as in amrclaw.

      integer max_mwaves
      parameter(max_mwaves = 10)

      integer method(7), mthlim(max_mwaves), mwaves, mcapa

      common /comclaw/ method, mthlim, mwaves, mcapa
