c     ==================
      subroutine setprob
c     ==================

      implicit double precision (a-h,o-z)
      common /comic/ qin(5),qout(5)
      common /param/  gamma,gamma1
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /cominf/ rinf,vinf,einf
c
c
c      # set idisc for cellave routines (see function fdisc)
       idisc = 2
c

c      # These should be read in as options
       gamma = 1.4d0
       gamma1 = gamma - 1.d0

       x0 = 0.5d0
       y0 = 0.0d0
       r0 = 0.2d0
       rhoin = 0.1d0
       pinf = 5.d0

c       read(7,*) gamma
c       gamma1 = gamma - 1.d0
c
cc      # read center and radius of bubble:
c       read(7,*) x0,y0,r0
cc      # density in bubble:
c       read(7,*) rhoin
cc      # pressure behind shock:
c       read(7,*) pinf
cc
c      # density outside bubble and pressure ahead of shock are fixed:
       rhoout = 1.d0
       pout   = 1.d0
       pin    = 1.d0

       qin(1) = rhoin
       qin(2) = 0.d0
       qin(3) = 0.d0
       qin(4) = pin/gamma1
       qin(5) = 1.d0

       qout(1) = rhoout
       qout(2) = 0.d0
       qout(3) = 0.d0
       qout(4) = pout/gamma1
       qout(5) = 0.d0
c
c     # Compute density and velocity behind shock from Hugoniot relations:
c     # ------------------------------------------------------------------

      rinf = ( gamma1 + (gamma+1)*pinf )/( (gamma+1) + gamma1*pinf )
      vinf = (1.0d0/sqrt(gamma)) * (pinf - 1.d0)/
     &       sqrt( 0.5*((gamma+1)/gamma) * pinf + 0.5*gamma1/gamma )
      einf = 0.5*rinf*vinf*vinf + pinf/gamma1

C       write(6,601) pinf,rinf,vinf,einf
C   601 format('pinf,rinf,vinf,einf:',/, 4e14.6)
C
      return
      end
