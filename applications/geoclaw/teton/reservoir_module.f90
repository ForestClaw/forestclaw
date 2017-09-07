MODULE reservoir_module
  IMPLICIT NONE
  INTEGER, PRIVATE, PARAMETER :: n = 32   !! number of points in the polygon
  DOUBLE COMPLEX, PRIVATE, DIMENSION(n) :: z, fe, dz    !! internal arrays
  DOUBLE PRECISION, PRIVATE, DIMENSION(2) :: xll,xur    !! corners of a box containing reservoir


CONTAINS
  LOGICAL FUNCTION point_in_poly(x0,y0)
    IMPLICIT NONE

    REAL(kind=8) :: x0,y0

    DOUBLE COMPLEX :: z0, trapsum
    REAL(kind=8) :: tol
    INTEGER :: i

    DOUBLE PRECISION :: pi

    pi = 4.d0*ATAN(1.d0)

    !! First, do simple check

    IF (x0 < xll(1) .OR. x0 > xur(1)) THEN
       point_in_poly = .FALSE.
       RETURN
    ENDIF

    IF (y0 < xll(2) .OR. y0 > xur(2)) THEN
       point_in_poly = .FALSE.
       RETURN
    ENDIF



    z0 = CMPLX(x0,y0,kind=8)



    !!n = 32

    DO i = 1,n
       !! z(i) = CMPLX(x(i),y(i),kind=8)
       fe(i) = 1.d0/(z(i) - z0)
!!       IF (i .GT. 1) THEN
!!          dz(i-1) = z(i) - z(i-1)
!!       ENDIF
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


  SUBROUTINE set_reservoir_path()
    IMPLICIT NONE

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: wpath
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wpath1d

    INTEGER :: i

    ALLOCATE(wpath(2,n))

    !! wpath1d automatically allocated
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
       z(i) = CMPLX(wpath(1,i),wpath(2,i),kind=8)
       IF (i .GT. 1) THEN
          dz(i-1) = z(i) - z(i-1)
       ENDIF
    END DO
    xll(1) = -111.6d0
    xur(1) = -111.2d0
    xll(2) = 43.8d0
    xur(2) = 44.0d0

    DEALLOCATE(wpath)

  END SUBROUTINE set_reservoir_path


end MODULE reservoir_module
