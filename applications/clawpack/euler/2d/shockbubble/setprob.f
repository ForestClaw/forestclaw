      subroutine shockbubble_setprob(gamma_in,
     &      x0_in,y0_in,r0_in,rhoin_in,pinf_in, idisc_in)

      implicit none
      double precision gamma_in, x0_in, y0_in, r0_in
      double precision rhoin_in, pinf_in
      double precision gamma,x0,y0,r0,rhoin,pinf
      double precision qin(5), qout(5)
      double precision gamma1
      double precision alf, beta
      double precision rinf,vinf,einf
      double precision rhoout, pout, pin
      integer idisc, idisc_in

      common /comic/ qin,qout
      common /cparam/  gamma,gamma1
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /cominf/ rinf,vinf,einf
c
c
c      # set idisc for cellave routines (see function fdisc)
       idisc = idisc_in
c

c      # These should be read in as options
       gamma = gamma_in
       gamma1 = gamma - 1.d0

       x0 = x0_in
       y0 = y0_in
       r0 = r0_in
       rhoin = rhoin_in
       pinf = pinf_in

c       x0 = 0.5d0
c       y0 = 0.0d0
c       r0 = 0.2d0
c       rhoin = 0.1d0
c       pinf = 5.d0

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
     &       sqrt( 0.5d0*((gamma+1)/gamma) * pinf + 0.5d0*gamma1/gamma )
      einf = 0.5d0*rinf*vinf*vinf + pinf/gamma1

      return
      end
