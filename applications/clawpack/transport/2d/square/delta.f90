double precision function delta(x)
    implicit none

    double precision x

    double precision pi, pi2
    common /compi/ pi, pi2

    double precision t0

    delta = exp(-x**2/(4*t0))/sqrt(4*pi*t0)

    return

end function delta

