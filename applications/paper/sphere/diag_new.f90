MODULE diag
  IMPLICIT NONE
  !
  ! Numerical mixing diagnostics are computed in Subroutine correlation_diag below.
  !
  ! Filament preservation diagnostic is computed in Subroutine filament_diag below.
  !
  ! It is assumed that the data is in one-dimensional arrays of length K, however,
  ! the code can easily be changed to acomodate other data structures.
  !
  ! dA(j) is the spherical area of grid cell j
  ! f1 is the mixing ratio of tracer 1
  ! f2 is the mixing ratio of tracer 2
  !

CONTAINS
  !
  ! compute correlation diagnostics
  !

  SUBROUTINE correlation_diag(f1,f2,K,dA,diagout)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: K
    REAL, DIMENSION(K)      , INTENT(IN) :: f1,f2,dA
    REAL :: sumA, diagout(3)
    !
    ! local workspace
    !
    REAL     :: root, tol,q1,q2,real_mixing,overshooting,range_pres_unmixing,c
    INTEGER  :: j
    REAL     :: q1_min,q1_max,q2_min,q2_max
    REAL     :: total_area, sqrt_arg

    q1_min = 0.1
    q1_max = 1.0

!! Changed so that end points of line segement are (q1_min,q2_min) and (q1_max,q2_max)
!! Note that "q2_min > q2_max".
    q2_min = corr_fct(q1_min)
    q2_max = corr_fct(q1_max)

    real_mixing          = 0.0
    overshooting         = 0.0
    range_pres_unmixing  = 0.0

    total_area = 0.0

    sumA = 0
    DO j = 1,K
       total_area = total_area+dA(j)

       q1 = f1(j)
       q2 = f2(j)

       !! Check to make sure we are not in the (very unlikely) situation where the argument
       !! to the sqrt (below) is negative.  This will happen if we are in a region which
       !! the cubic has three real roots, and so we'd have to be careful about which one we pick.
       sqrt_arg = (-DBLE(1687296) + DBLE(12168000)*q2 &
            -DBLE(29250000)*q2**2+DBLE(23437500)*q2**3+DBLE(29648025)*q1**2)
       IF (sqrt_arg < 0) THEN
          WRITE(6,*) 'Warning : (xk,yk) data is in region where there are three possible ', &
               ' real roots to the closest point problem'
          STOP
       ENDIF

       c = (DBLE(65340)*q1+12.*SQRT(-DBLE(1687296)+DBLE(12168000)*q2 &
            -DBLE(29250000)*q2**2+DBLE(23437500)*q2**3+DBLE(29648025)*q1**2))**(1./3.)
       c=c/(DBLE(60))

       root = c-(-(DBLE(13)/DBLE(75))+(DBLE(5)/DBLE(12))*q2)/c
       root = MAX(0.1,root)
       root = MIN(1.0,root)

       !! Also fixed bug here : call to line_fct passed in q2 instead of q1 as argument.
       IF (q2 < corr_fct(q1).AND.q2 > line_fct(q1,q1_min,q1_max,q2_min,q2_max)) THEN
          !
          ! `real' mixing
          !
          real_mixing = real_mixing + dist_fct(root,q1,q2)*dA(j)
       ELSE IF (q1 < q1_max.AND.q1 > q1_min.AND.q2 < q2_min.AND.q2 > q2_max) THEN
          !! Note that "q2_min > q2_max", so this 'if' branch had to be modifed

          !
          ! range-preserving unmixing
          !
          range_pres_unmixing = range_pres_unmixing+dist_fct(root,q1,q2)*dA(j)
       ELSE
          !
          ! overshooting
          !
          overshooting = overshooting + dist_fct(root,q1,q2)*dA(j)
       END IF
       sumA = sumA + dist_fct(root,q1,q2)*dA(j)
    END DO
    WRITE(*,*) "========================================================================"
    WRITE(*,*) " "
    WRITE(*,*) "mixing diagnostics"
    WRITE(*,*) " "
    WRITE(*,*) "------------------------------------------------------------------------"
    WRITE(*,*) " "
    WRITE(*,*) "real_mixing ",real_mixing/total_area
    WRITE(*,*) "range_pres_unmixing ",range_pres_unmixing/total_area
    WRITE(*,*) "overshooting     ",overshooting/total_area
    WRITE(*,*) " "
    WRITE(*,*) "real_mixing (frac) ", (real_mixing/sumA)
    WRITE(*,*) "range_pres_mixing (frac) ", (range_pres_unmixing/sumA)
    WRITE(*,*) "over_shooting (frac) ", (overshooting/sumA)
    WRITE(*,*) "sumA " , sumA
    WRITE(*,*) "========================================================================"

    diagout(1) = real_mixing/total_area
    diagout(2) = range_pres_unmixing/total_area
    diagout(3) = overshooting/total_area
  END SUBROUTINE correlation_diag

  !
  ! correlation function
  !
  REAL FUNCTION corr_fct(x)
    IMPLICIT NONE
    REAL , INTENT(IN)  :: x
    corr_fct = -0.8*x**2+0.9
  END FUNCTION corr_fct
  !
  ! Eucledian distance function
  !
  REAL FUNCTION dist_fct(x,x0,y0)
    IMPLICIT NONE
    REAL , INTENT(IN)  :: x,x0,y0
    dist_fct = SQRT((x-x0)*(x-x0)/(0.9**2)+(corr_fct(x)-y0)*(corr_fct(x)-y0)/(0.792**2))
  END FUNCTION dist_fct
  !
  ! straight line line function
  !
  REAL FUNCTION line_fct(x,xmin,xmax,ymin,ymax)
    IMPLICIT NONE
    REAL , INTENT(IN)  :: x,xmin,xmax,ymin,ymax
    REAL  :: a,b
    !
    ! line: y=a*x+b
    !
    a = (ymax-ymin)/(xmax-xmin)
    b = ymin-xmin*a
    line_fct = a*x+b
  END FUNCTION line_fct

  !
  ! linit = .TRUE. if t=0 else .FALSE.
  !

SUBROUTINE filament_diag(K,f1,dA,fila_t0,linit)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                    :: K
    REAL(kind=8)   , DIMENSION(K)  , INTENT(IN)    :: f1, dA
    REAL(kind=8)   , DIMENSION(100), INTENT(INOUT) :: fila_t0
    LOGICAL        , INTENT(IN)    :: linit
    !
    ! local workspace
    !
    REAL(kind=8) :: threshold,tiny
    REAL(kind=8) :: area_out, a, area_in
    INTEGER :: j,jk,jlevels

    double precision b_init

    tiny = 1.0d-8
    b_init = 0.0d0
    
    OPEN (unit = 31, file='filament.dat',status='replace')
    jlevels = 20
    DO jk = 0,jlevels
       threshold = b_init + jk/20.d0
       area_out = 0.0
       area_in = 0
       DO j = 1,K
          IF (f1(j) .GE. threshold-tiny) THEN
             area_out = area_out + dA(j)
          else
             area_in = area_in + dA(j)
          END IF
       END DO
       IF (linit) THEN
          fila_t0(jk) = area_out
       ELSE
          IF (fila_t0(jk) .EQ. 0) THEN
             a = 0.d0
          ELSE
             a = 100*area_out/fila_t0(jk)
          ENDIF
          WRITE(31,100) threshold, a, fila_t0(jk)
       END IF
       write(6,110) threshold,area_out, area_in, area_out+area_in, fila_t0(jk)
    END DO
    CLOSE(31)
100 FORMAT(F12.2,F12.4, F24.16)
110 FORMAT(F12.2,5F24.16)
  END SUBROUTINE filament_diag
END MODULE diag
