      subroutine setprob()
      implicit none

      double precision rho, bulk, cc, zz
      common /cparam/ rho,bulk,cc,zz

      double precision pi, pi2
      common /compi/ pi, pi2

      pi = 4.d0*atan(1.d0)
      pi2 = 2*pi

c     # These need to be assigned from user options
      open(10,file='setprob.data')
      read(10,*) rho
      read(10,*) bulk
      close(10)

      cc = sqrt(bulk/rho)
      zz = rho*cc

      return
      end
