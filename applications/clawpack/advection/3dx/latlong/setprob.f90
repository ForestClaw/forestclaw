subroutine setprob()
    implicit none

    double precision pi, pi2
    common /compi/ pi, pi2

    integer example
    common /common_ex/ example

    double precision revs_per_second
    common /com_latlong/ revs_per_second

    integer manifold
    common /com_manifold/ manifold

    double precision maxelev
    common /com_extruded/ maxelev

    double precision longitude(2), latitude(2)

    pi = 4.d0*atan(1.d0)
    pi2 = 2.0*pi

    open(10,file='setprob.data')
    read(10,*) example
    read(10,*) manifold
    read(10,*) revs_per_second
    read(10,*) longitude(1)
    read(10,*) longitude(2)
    read(10,*) latitude(1)
    read(10,*) latitude(2)
    read(10,*) maxelev
    close(10)

end
