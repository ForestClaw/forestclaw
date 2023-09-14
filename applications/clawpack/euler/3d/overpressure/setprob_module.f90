module setprob_mod
    implicit none

    double precision pi, pi2
    integer example
    integer mapping
    integer manifold, mcapa
    integer init_choice
    double precision gamma, gamma1
    double precision rhoin, rhoout, pin, pout
    double precision x0, y0, z0, r0, maxelev
    double precision qin(5), qout(5)

end module setprob_mod

subroutine setprob
    use setprob_mod

    implicit none

    !! Not needed in a common blocks
    double precision longitude(2), latitude(2)

    double precision pi_com, pi2_com
    common /compi/ pi_com, pi2_com

    integer mapping_com
    common /com_mapping/ mapping_com


    pi = 4.d0*atan(1.d0)
    pi2 = 2*pi

    pi_com = pi
    pi2_com = pi2

    open(10,file='setprob.data')
    read(10,*) example
    read(10,*) mapping
    read(10,*) manifold
    read(10,*) mcapa
    read(10,*) init_choice

    mapping_com = mapping

    !! # These should be read in as options
    read(10,*) gamma
    gamma1 = gamma - 1.d0

    read(10,*) x0    
    read(10,*) y0    
    read(10,*) z0
    read(10,*) r0    
    read(10,*) rhoin 
    read(10,*) rhoout
    read(10,*) pin
    read(10,*) pout
    read(10,*) longitude(1)
    read(10,*) longitude(2)
    read(10,*) latitude(1)
    read(10,*) latitude(2)
    read(10,*) maxelev
    close(10)


    !! # density outside bubble and pressure ahead of shock are fixed:

    qin(1) = rhoin
    qin(2) = 0.d0
    qin(3) = 0.d0
    qin(4) = 0.d0
    qin(5) = pin/gamma1

    qout(1) = rhoout
    qout(2) = 0.d0
    qout(3) = 0.d0
    qout(4) = 0.d0
    qout(5) = pout/gamma1

    return
end subroutine setprob
