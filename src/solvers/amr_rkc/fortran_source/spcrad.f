      double precision function spcrad(neqn,t,y)
      implicit none
      integer          neqn
      double precision t,y(neqn)

c     # Estimate the spectral radius of the operator.
      spcrad = 0

      return
      end
