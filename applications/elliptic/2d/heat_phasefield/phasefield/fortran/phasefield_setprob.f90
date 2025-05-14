SUBROUTINE phasefield_setprob()
    IMPLICIT none

    INTEGER example
    COMMON /comm_example/ example

    DOUBLE PRECISION S, alpha, m, xi, k, gamma
    COMMON /comm_parms/ S, alpha, m, xi, k, gamma

    DOUBLE PRECISION r0, x0, y0
    COMMON /comm_init/ r0, x0, y0

    INTEGER bc_type(0:3)
    COMMON /comm_bc/ bc_type

    DOUBLE PRECISION pi,pi2
    COMMON /compi/ pi, pi2

    pi = 4.d0*atan(1.d0)
    pi2 = 2*pi

    open(10,file='setprob.data')
    read(10,*) example
    read(10,*) S
    read(10,*) alpha
    read(10,*) m
    read(10,*) xi
    read(10,*) k
    read(10,*) gamma
    read(10,*) r0
    read(10,*) x0
    read(10,*) y0

    read(10,*) bc_type(0)
    read(10,*) bc_type(1)
    read(10,*) bc_type(2)
    read(10,*) bc_type(3)
    close(10)

end subroutine phasefield_setprob