SUBROUTINE setprob()
    IMPLICIT NONE

    INTEGER iunit
    CHARACTER(len=25) fname      

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION :: grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level

    DOUBLE PRECISION :: alpha, breaking
    COMMON /common_dispersive/ alpha, breaking

    DOUBLE PRECISION :: a,b,h0
    COMMON /common_initheight/ a,b,h0

    pi = 4.d0*atan(1.d0)
    pi2 = 2*pi

    iunit = 10
    fname = 'setprob.data'      
    call opendatafile(iunit, fname)
    read(10,*) grav
    read(10,*) dry_tolerance
    read(10,*) sea_level

    !! Lannes method
    read(10,*) breaking
    read(10,*) alpha

    !! Initial conditions
    read(10,*) a
    read(10,*) b
    read(10,*) h0
    close(10);

end subroutine setprob
