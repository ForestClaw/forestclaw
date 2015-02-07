      subroutine setprob()
      implicit none
      double precision rho, bulk, cc, zz
      common /cparam/ rho,bulk,cc,zz

c     # These need to be assigned from user options
      rho = 1.d0
      bulk = 4.0
      cc = sqrt(bulk/rho)
      zz = rho*cc

      return
      end
