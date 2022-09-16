subroutine setprob()
    implicit none

    double precision pi, pi2
    common /compi/ pi, pi2

    integer example
    common /common_ex/ example

    integer manifold
    common /com_manifold/ manifold

    double precision revs_per_second
    common /com_latlong/ revs_per_second

    pi = 4.d0*atan(1.d0)
    pi2 = 2.0*pi

    open(10,file='setprob.data')
    read(10,*) example
    read(10,*) manifold
    read(10,*) revs_per_second
    !! latitude, longitude
    close(10)



end
