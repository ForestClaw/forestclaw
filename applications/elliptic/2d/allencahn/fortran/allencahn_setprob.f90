subroutine allencahn_setprob()
    use hsmooth_mod
    implicit none

    INTEGER example
    COMMON /comm_example/ example

    INTEGER beta_choice
    COMMON /comm_beta/ beta_choice

    DOUBLE PRECISION alpha, x0,y0,a,b
    COMMON /comm_rhs/ alpha, x0,y0,a,b

    integer bc_type(0:3)
    common /comm_bc/ bc_type

    DOUBLE PRECISION pi,pi2
    COMMON /compi/ pi, pi2

    integer i

    pi = 4.d0*atan(1.d0)
    pi2 = 2*pi

    open(10,file='setprob.data')
    read(10,*) example
    read(10,*) beta_choice
    read(10,*) alpha
    read(10,*) x0
    read(10,*) y0
    read(10,*) a
    read(10,*) b
    read(10,*) eps_disk
    read(10,*) m_polar

    call allocate_polar_arrays()
    read(10,*) (x0_polar(i),i=1,m_polar)
    read(10,*) (y0_polar(i),i=1,m_polar)
    read(10,*) (r0_polar(i),i=1,m_polar)
    read(10,*) (r1_polar(i),i=1,m_polar)
    read(10,*) (n_polar(i) ,i=1,m_polar)

    read(10,*) bc_type(0)
    read(10,*) bc_type(1)
    read(10,*) bc_type(2)
    read(10,*) bc_type(3)
    close(10)

end subroutine allencahn_setprob