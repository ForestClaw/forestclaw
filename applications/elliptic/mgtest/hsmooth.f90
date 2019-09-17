module hsmooth_mod
    implicit none
    integer :: m_polar = -1

    double precision, dimension(:), allocatable :: x0_polar,y0_polar
    double precision, dimension(:), allocatable :: r0_polar,r1_polar
    integer, dimension(:), allocatable :: n_polar

    double precision eps_disk

contains
    subroutine allocate_polar_arrays()
        implicit none

        if (m_polar .lt. 0) then
            write(6,*) 'hsmooth:allocate_polar_array : m_polar is not defined'
            stop
        endif

        allocate(x0_polar(m_polar), y0_polar(m_polar),r0_polar(m_polar), & 
                 r1_polar(m_polar), n_polar(m_polar))
    end subroutine allocate_polar_arrays
end module hsmooth_mod


double precision function polar_interface(id,theta)
    implicit none

    double precision theta
    integer id

    double precision p, dpdt, d2pdt2

    call polar_interface_complete(id,theta, p, dPdt, d2pdt2)
    polar_interface = p

end function polar_interface

subroutine polar_interface_complete(id,theta, p, dPdtheta, d2Pdtheta2)
    use hsmooth_mod, only : r0_polar, r1_polar, n_polar
    implicit none

    double precision theta, P, dPdtheta, d2Pdtheta2
    integer id

    double precision r0, r1, n

    r0 = r0_polar(id)
    r1 = r1_polar(id)
    n = n_polar(id)    

    p = r0*(1 + r1*cos(n*theta))
    dpdtheta   = r0*(-n*r1*sin(n*theta))
    d2Pdtheta2 = r0*(-n**2*r1*cos(n*theta))

end subroutine polar_interface_complete

double precision function Hsmooth(id,r,theta)
    use hsmooth_mod, only : eps_disk
    implicit none

    double precision r, theta
    integer id

    double precision p, polar_interface, a

    a = eps_disk

    P = polar_interface(id,theta)
    Hsmooth = (tanh((r-P)/a) + 1)/2

end function Hsmooth

subroutine Hsmooth_grad(id,r,theta,grad_H)
    use hsmooth_mod, only : eps_disk
    implicit none

    double precision r, theta, grad_H(2)
    integer id

    double precision a, a2, sech, sech2, p, dpdt, dp2dt2

    if (r .eq. 0) then
        r = 1d-15
    endif

    a = eps_disk
    a2 = 2*a

    call polar_interface_complete(id,theta,P,dpdt,dp2dt2)

    sech2 = sech((r-P)/a)**2
    grad_H(1) = sech2/a2
    grad_H(2) = -dpdt*sech2/(a2*r)

end subroutine Hsmooth_grad

double precision function Hsmooth_laplacian(id,r,theta)
    use hsmooth_mod, only : eps_disk
    implicit none

    double precision r, theta
    integer id

    double precision sech, a, a2, r2, p, dpdt, d2pdt2
    double precision s1, s2, s3, s4, t, sech2, st

    if (r .eq. 0) then
        r = 1d-15
    endif

    a = eps_disk
    a2 = a**2
    r2 = r**2    

    call polar_interface_complete(id,theta,p,dpdt,d2pdt2)

    sech2 = sech((r-p)/a)**2
    t = tanh((r-p)/a)
    st = t*sech2
    s1 = dpdt**2*st/a2
    s2 = d2pdt2*sech2/(2*a)
    s3 = st/a2
    s4 = sech2/(2*a*r)
    Hsmooth_laplacian = (-s1-s2)/r2 - s3 + s4

end function Hsmooth_laplacian


