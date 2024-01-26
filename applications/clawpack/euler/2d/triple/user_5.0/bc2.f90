SUBROUTINE clawpack5_bc2(meqn,mbc,mx,my,xlower,ylower, &
     dx,dy,q,maux,aux,t,dt,mthbc)

  IMPLICIT NONE
  INTEGER meqn,mbc,mx,my,maux,mthbc(4)
  DOUBLE PRECISION xlower,ylower,dx,dy,t,dt

  DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

  INTEGER i,j,ibc,jbc,m


  !!-------------------------------------------------------
  !!     # left boundary:
  !!-------------------------------------------------------
  GOTO (100,110,120,130) mthbc(1)+1
!! this is how we skip over this side... if (mthbc(1)+1
!! is not 1,2,3 or 4, then the goto above falls through to here...
  goto 199

  100 continue
!! # user-specified boundary conditions go here in place of error output
  WRITE(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc2'
  STOP
  go to 199

  110 continue
  !! # zero-order extrapolation:
  DO j = 1-mbc, my+mbc
     DO ibc=1,mbc
        DO m=1,meqn
           q(m,1-ibc,j) = q(m,1,j)
        ENDDO
     ENDDO
  ENDDO
  go to 199

  120 continue
  !! # periodic:
  DO j = 1-mbc, my+mbc
     DO ibc=1,mbc
        DO m=1,meqn
           q(m,1-ibc,j) = q(m,mx+1-ibc,j)
        ENDDO
     ENDDO
  ENDDO
  GOTO 199

  130 continue
  !! # solid wall (assumes 2'nd component is velocity or momentum in x):
  DO j = 1-mbc, my+mbc
     DO ibc=1,mbc
        DO m=1,meqn
           q(m,1-ibc,j) = q(m,ibc,j)
        ENDDO
     ENDDO
  ENDDO

  !! # negate the normal velocity:
  DO j = 1-mbc, my+mbc
     DO ibc=1,mbc
        q(2,1-ibc,j) = -q(2,ibc,j)
     ENDDO
  ENDDO
  go to 199

  199 continue

  !!-------------------------------------------------------
  !!     # right boundary:
  !!-------------------------------------------------------
  GOTO (200,210,220,230) mthbc(2)+1
  GOTO 299

200 CONTINUE
  !! # user-specified boundary conditions go here in place of error output
  WRITE(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
  STOP
  go to 299

  210 continue
  !! # zero-order extrapolation:
  DO j = 1-mbc, my+mbc
     DO ibc=1,mbc
        DO m=1,meqn
           q(m,mx+ibc,j) = q(m,mx,j)
        ENDDO
     ENDDO
  ENDDO
  go to 299

220 CONTINUE
  !! # periodic:
  DO j = 1-mbc, my+mbc
     DO ibc=1,mbc
        DO m=1,meqn
           q(m,mx+ibc,j) = q(m,ibc,j)
        ENDDO
     ENDDO
  ENDDO
  go to 299

  230 continue
  !! # solid wall (assumes 2'nd component is velocity or momentum in x):
  DO j = 1-mbc, my+mbc
     DO ibc=1,mbc
        DO m=1,meqn
           q(m,mx+ibc,j) = q(m,mx+1-ibc,j)
        ENDDO
     ENDDO
  ENDDO

  !! # negate the normal velocity:
  DO j = 1-mbc, my+mbc
     DO ibc=1,mbc
        q(2,mx+ibc,j) = -q(2,mx+1-ibc,j)
     ENDDO
  ENDDO
  go to 299

299 CONTINUE

  !!-------------------------------------------------------
  !!     # bottom boundary:
  !!-------------------------------------------------------
  GOTO (300,310,320,330) mthbc(3)+1
  GOTO 399

300 CONTINUE
  !! # user-specified boundary conditions go here in place of error output
  WRITE(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
  STOP
  go to 399

310 CONTINUE
  !! # zero-order extrapolation:
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        DO m=1,meqn
           q(m,i,1-jbc) = q(m,i,1)
        ENDDO
     ENDDO
  ENDDO
  go to 399

320 CONTINUE
  !! # periodic:
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        DO m=1,meqn
           q(m,i,1-jbc) = q(m,i,my+1-jbc)
        ENDDO
     ENDDO
  ENDDO
  go to 399

330 CONTINUE
  !! # solid wall (assumes 3'rd component is velocity or momentum in y):
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        DO m=1,meqn
           q(m,i,1-jbc) = q(m,i,jbc)
        ENDDO
     ENDDO
  ENDDO

  !! # negate the normal velocity:
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        q(3,i,1-jbc) = -q(3,i,jbc)
     ENDDO
  ENDDO
  go to 399

399 CONTINUE

  !!-------------------------------------------------------
  !!     # top boundary:
  !!-------------------------------------------------------
  GOTO (400,410,420,430) mthbc(4)+1
  GOTO 499

400 CONTINUE
!! # user-specified boundary conditions go here in place of error output
  WRITE(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
  STOP
  go to 499

410 CONTINUE
  !! # zero-order extrapolation:
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        DO m=1,meqn
           q(m,i,my+jbc) = q(m,i,my)
        ENDDO
     ENDDO
  ENDDO
  go to 499

420 CONTINUE
  !! # periodic:
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        DO m=1,meqn
           q(m,i,my+jbc) = q(m,i,jbc)
        ENDDO
     ENDDO
  ENDDO
  go to 499

  430 continue
  !! # solid wall (assumes 3'rd component is velocity or momentum in y):
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        DO m=1,meqn
           q(m,i,my+jbc) = q(m,i,my+1-jbc)
        ENDDO
     ENDDO
  ENDDO

  !! # negate the normal velocity:
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        q(3,i,my+jbc) = -q(3,i,my+1-jbc)
     ENDDO
  ENDDO
  go to 499

499 CONTINUE

  RETURN
END SUBROUTINE clawpack5_bc2
