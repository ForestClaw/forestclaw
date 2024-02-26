double precision function  fdisc(blockno,xc,yc,zc)
    implicit none

    double precision :: xc,yc,zc, xp, yp, zp
    integer :: blockno
    integer*8 :: cont, fclaw_map_get_context

    double precision :: r

    cont = fclaw_map_get_context()

    call fclaw_map_3d_c2m(cont,blockno,xc,yc,zc,xp,yp,zp)

    r = sqrt((xp-0.5d0)**2 + (yp-1.d0)**2)

    fdisc = r-0.25d0
end function fdisc
