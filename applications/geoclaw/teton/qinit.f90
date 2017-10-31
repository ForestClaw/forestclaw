! qinit routine for parabolic bowl problem, only single layer
SUBROUTINE teton_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use geoclaw_module, only: grav
    use reservoir_module

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Parameters for problem
    REAL(kind=8), PARAMETER :: a = 1.d0
    REAL(kind=8), PARAMETER :: sigma = 0.5d0

    !! Original code`
    !! REAL(kind=8), PARAMETER :: xll = -112.3895d0
    !! REAL(kind=8), PARAMETER :: x0 = -111.5391666666667d0  !! x-coords of Teton Dam
    !! REAL(kind=8), PARAMETER :: x1 = -111.24d0  !! Left edge of domain
    !! REAL(kind=8), PARAMETER :: h0 = 1540.d0
    !! REAL(kind=8), PARAMETER :: h1 = 1720.d0


    !! 9/4/2017 : Using larger domain
    REAL(kind=8), PARAMETER :: xll = -112.34626736d0
    REAL(kind=8), PARAMETER :: x0 = -111.5391666666667d0  !! x-coords of Teton Dam
    REAL(kind=8), PARAMETER :: x1 = -111.24d0  !! Left edge of domain
    REAL(kind=8), PARAMETER :: h0 = 1560.d0
    REAL(kind=8), PARAMETER :: h1 = 1720.d0


    ! Other storage
    INTEGER :: i,j
    REAL(kind=8) :: omega,x,y,eta

    !! Capacity of reservoir at time of release : 251000 acre-feet
    !! 1 acre-foot = 1233.4829 m^3
    !! Conversion of deg2m2 = 8948239546.22

    !! Actual mass : 251000*1233.4829 = 309,604,207 (m^3)

    !! Total mass as reported by GeoClaw (using h0=1615, h1=1700)
    !!        is  27,298,685 (m^3)  (assuming mass reported by GeoClaw is in m^3)
    !! Total mass as reported by GeoClaw (using h0=1625.d0, h1 = 1700)
    !!        is 62,069,921  (m^3)
    !! 9/4/2017 (using deg2m2 conversion)
    !! Total mass as reported by GeoClaw (using h0=1540, h1 = 1720, maxlevel=6)
    !!      is mass = 263,167,751.82 m^3
    !! Total mass as reported by GeoClaw (using h0=1540, h1 = 1720, maxlevel=7)
    !!      is mass = 274,799,319.07 m^3
 
    CALL set_reservoir_path()

    DO i=1-mbc,mx+mbc
        x = xlower + (i - 0.5d0)*dx
        do j=1-mbc,my+mbc
           y = ylower + (j - 0.5d0) * dx
           q(1,i,j) = 0
           IF (x .GT. x0) THEN
              IF (point_in_poly(x,y)) THEN
                 eta = h0 + (h1-h0)/(x1-x0)*(x - x0)
                 q(1,i,j) =  MAX(0.d0,eta - aux(1,i,j))
              ENDIF
           ENDIF
           q(2,i,j) = 0.d0
           q(3,i,j) = 0.d0
        enddo
     ENDDO

   END SUBROUTINE teton_qinit


