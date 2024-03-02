double precision function fdisc(blockno,xc,yc)
    implicit none

    integer :: blockno
    double precision :: xc,yc

    double precision :: x0, y0, r0
    common/cdisc/ x0,y0, r0


    integer*8 :: cont, get_context
    double precision :: f1,f2

    !! # circle of radius r0:
    fdisc = (xc-x0)**2 + (yc-y0)**2 - r0**2
    return

end function fdisc
