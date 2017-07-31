! qinit routine for parabolic bowl problem, only single layer
SUBROUTINE teton_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use geoclaw_module, only: grav

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Parameters for problem
    REAL(kind=8), PARAMETER :: a = 1.d0
    REAL(kind=8), PARAMETER :: sigma = 0.5d0

    REAL(kind=8), PARAMETER :: xll = -112.3895d0
    REAL(kind=8), PARAMETER :: x0 = -111.5391666666667d0  !! x-coords of Teton Dam
    REAL(kind=8), PARAMETER :: x1 = -111.24d0  !! Left edge of domain
    REAL(kind=8), PARAMETER :: h0 = 1540.d0
    REAL(kind=8), PARAMETER :: h1 = 1720.d0
!!    REAL(kind=8), PARAMETER :: h0 = 1625.d0   !! Flat surface
!!    REAL(kind=8), PARAMETER :: h1 = 1720.d0

    INTEGER, PARAMETER :: nmax = 100
    REAL(kind=8), DIMENSION(nmax) :: xp,yp
    LOGICAL point_in_poly

    ! Other storage
    INTEGER :: i,j
    REAL(kind=8) :: omega,x,y,eta

    !! Capacity of reservoir at time of release : 251000 acre-feet
    !! 1 acre-foot = 1233.4829 m^3
    !! Actual mass : 251000*1233.4829 = 309,604,207 (m^3)
    !! Total mass as reported by GeoClaw (using h0=1615, h1=1700)
    !!        is  27,298,685 (m^3)  (assuming mass reported by GeoClaw is in m^3)
    !! Total mass as reported by GeoClaw (using h0=1625.d0, h1 = 1700)
    !!        is 62,069,921  (m^3)

    CALL set_reservoir_path(nmax,xp,yp)

    DO i=1-mbc,mx+mbc
        x = xlower + (i - 0.5d0)*dx
        do j=1-mbc,my+mbc
           y = ylower + (j - 0.5d0) * dx
           q(1,i,j) = 0
           IF (x .GT. x0) THEN
              IF (point_in_poly(nmax,xp,yp,x,y)) THEN
                 eta = h0 + (h1-h0)/(x1-x0)*(x - x0)
                 q(1,i,j) =  MAX(0.d0,eta - aux(1,i,j))
              ENDIF
           ENDIF
           q(2,i,j) = 0.d0
           q(3,i,j) = 0.d0
        enddo
     ENDDO

   END SUBROUTINE teton_qinit


  LOGICAL FUNCTION point_in_poly(nmax,x,y,x0,y0)
    IMPLICIT NONE

    INTEGER :: nmax
    REAL(kind=8), DIMENSION(nmax) :: x,y
    REAL(kind=8) :: x0,y0

    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: z, fe, dz
    DOUBLE COMPLEX :: z0, trapsum
    REAL(kind=8) :: tol
    INTEGER :: i,n

    DOUBLE PRECISION :: pi

    pi = 4.d0*ATAN(1.d0)

    z0 = CMPLX(x0,y0,kind=8)

    n = 32
    ALLOCATE(z(n), fe(n), dz(n))

    DO i = 1,n
       z(i) = CMPLX(x(i),y(i),kind=8)
       fe(i) = 1.d0/(z(i) - z0)
       IF (i .GT. 1) THEN
          dz(i-1) = z(i) - z(i-1)
       ENDIF
    ENDDO
100 FORMAT(I3, 4E24.16)

    !! Use trapezoidal rule to compute winding number
    trapsum = CMPLX(0.d0,0.d0,kind=8)
    DO i = 1,n-1
       trapsum = trapsum + (fe(i) + fe(i+1))*dz(i)/2.d0
    ENDDO
    trapsum = trapsum/(2*pi*CMPLX(0,1.d0,kind=8))

    tol = 0.25d0
    IF (ABS(trapsum) .LT. tol) THEN
       point_in_poly = .FALSE.
    ELSE
       point_in_poly = .TRUE.
    ENDIF


  END FUNCTION point_in_poly


SUBROUTINE set_reservoir_path(nmax,xpath,ypath)
  IMPLICIT NONE

  INTEGER nmax
  DOUBLE PRECISION, DIMENSION(nmax) :: xpath,ypath
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: wpath
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wpath1d

  INTEGER :: n, i

  n = 32
  ALLOCATE(wpath(2,n))

  wpath1d = (/-111.5329602010579d0,43.90768251711403d0,  &
       -111.5425694320295d0,  43.91291520283087d0,  &
       -111.5265517148000d0,  43.92985173935236d0,     &
       -111.4838301614854d0,  43.93342747699161d0,  &
       -111.4526923402719d0,  43.92275343741503d0,  &
       -111.4229733394132d0,  43.92102950267466d0,  &
       -111.3888345969115d0,  43.94256354582464d0,  &
       -111.3547747442403d0,  43.96012412188881d0,  &
       -111.3290981366292d0,  43.94067321304157d0,  &
       -111.3055308323582d0,  43.92361145688078d0,   &
       -111.2787647823537d0,  43.93643579897589d0,   &
       -111.2483464534973d0,  43.94374241316081d0,   &
       -111.2512048905199d0,  43.93181718339741d0,   &
       -111.3004354725815d0,  43.91434621684393d0,   &
       -111.3310923535663d0,  43.92063116805807d0,   &
       -111.3574777595737d0,  43.93748955926343d0,   &
       -111.3734430072187d0,  43.92390465892466d0,   &
       -111.3830291449872d0,  43.91639601995375d0,   &
       -111.4070584125505d0,  43.91053159674234d0,   &
       -111.4273298964650d0,  43.90372138216379d0,    &
       -111.4497771036826d0,  43.90847035302123d0,   &
       -111.4470054376565d0,  43.89319391079540d0,    &
       -111.4340292718339d0,  43.88128609194279d0,   &
       -111.4328310870467d0,  43.87453080137639d0,   &
       -111.4525518180650d0,  43.87284083784596d0,    &
       -111.4655146366021d0,  43.88648791247213d0,   &
       -111.4742548002720d0,  43.89924049270922d0,    &
       -111.4743486255008d0,  43.91523470848247d0,   &
       -111.4878354873365d0,  43.92113549058154d0,   &
       -111.5039146384135d0,  43.91948128966648d0,   &
       -111.5182133156775d0,  43.91450172601051d0,   &
       -111.5329602010579d0,  43.90768251711403d0/)   !! Last line repeats first line

  wpath = RESHAPE(wpath1d,(/2,n/))

  DO i = 1,n
     xpath(i) = wpath(1,i)
     ypath(i) = wpath(2,i)
  END DO

END SUBROUTINE set_reservoir_path
