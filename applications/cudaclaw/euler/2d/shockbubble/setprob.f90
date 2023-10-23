subroutine setprob()
    implicit none

    double precision :: qin(5), qout(5)
    common /comic/ qin,qout

    double precision :: gamma, gamma1
    common /cparam/  gamma,gamma1

    integer :: idisc
    double precision :: x0, y0, alf, beta, r0
    common/cdisc/ x0,y0,alf,beta,r0,idisc

    double precision :: rinf, vinf, einf
    common /cominf/ rinf,vinf,einf

    double precision :: rhoin, rhoout, pout, pin, pinf

    open(10,file='setprob.data')
       
    !! # These should be read in as options
    read(10,*) gamma
    gamma1 = gamma - 1.d0

    read(10,*) x0    
    read(10,*) y0    
    read(10,*) r0    
    read(10,*) rhoin 
    read(10,*) pinf  

    !! # set idisc for cellave routines (see function fdisc)
    read(10,*) idisc

    close(10)


    !! # density outside bubble and pressure ahead of shock are fixed:
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

    !! # Compute density and velocity behind shock from Hugoniot relations:
    !! # ------------------------------------------------------------------

    rinf = ( gamma1 + (gamma+1)*pinf )/( (gamma+1) + gamma1*pinf )
    vinf = (1.0d0/sqrt(gamma)) * (pinf - 1.d0)/ & 
             sqrt( 0.5d0*((gamma+1)/gamma) * pinf + 0.5d0*gamma1/gamma )
    einf = 0.5d0*rinf*vinf*vinf + pinf/gamma1

    return
end subroutine setprob
