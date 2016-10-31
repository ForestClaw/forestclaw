SUBROUTINE clawapck5_bc2(meqn,mbc,mx,my,xlower,ylower, &
     dx,dy,q,maux,aux,t,dt,mthbc)
  !!     =====================================================
  !!
  !!     # Standard boundary condition choices for claw2
  !!
  !!     # At each boundary  k = 1 (left),  2 (right),  3 (top), 4 (bottom):
  !!     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
  !!     #            =  1  for zero-order extrapolation
  !!     #            =  2  for periodic boundary coniditions
  !!     #            =  3  for solid walls, assuming this can be implemented
  !!     #                  by reflecting the data about the boundary and then
  !!     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
  !!     #                  component of q.
  !!     ------------------------------------------------
  !!
  !!     # Extend the data from the interior cells (1:mx, 1:my)
  !!     # to the ghost cells outside the region:
  !!     #   (i, 1-jbc)   for jbc = 1,mbc,  i = 1-mbc, mx+mbc
  !!     #   (i, my+jbc)  for jbc = 1,mbc,  i = 1-mbc, mx+mbc
  !!     #   (1-ibc, j)   for ibc = 1,mbc,  j = 1-mbc, my+mbc
  !!     #   (mx+ibc, j)  for ibc = 1,mbc,  j = 1-mbc, my+mbc
  !!
  implicit double precision (a-h,o-z)
  dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
  dimension mthbc(4)

  !!
  !!
  !!-------------------------------------------------------
  !!     # left boundary:
  !!-------------------------------------------------------
  go to (100,110,120,130) mthbc(1)+1

100 continue
  !!     # user-specified boundary conditions go here in place of error output
  write(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc2'
  stop
  go to 199

110 continue
  !! # zero-order extrapolation:
  DO j = 1-mbc, my+mbc
     DO ibc=1,mbc
        DO m=1,meqn
           q(m,1-ibc,j) = q(m,1,j)
115     ENDDO
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
  go to 199

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
  DO 136 j = 1-mbc, my+mbc
     DO 136 ibc=1,mbc
        q(2,1-ibc,j) = -q(2,ibc,j)
     ENDDO
  ENDDO
  go to 199

199 continue
  !!
  !!-------------------------------------------------------
  !!     # right boundary:
  !!-------------------------------------------------------
  go to (200,210,220,230) mthbc(2)+1

200 continue
  !! # user-specified boundary conditions go here in place of error output
  write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
  stop
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

220 continue
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

299 continue
  !!
  !!-------------------------------------------------------
  !!     # bottom boundary:
  !!-------------------------------------------------------
  go to (300,310,320,330) mthbc(3)+1

300 CONTINUE
  !! # user-specified boundary conditions go here in place of error output
  write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
  stop
  go to 399

310 continue
  !! # zero-order extrapolation:
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        DO m=1,meqn
           q(m,i,1-jbc) = q(m,i,1)
        ENDDO
     ENDDO
  ENDDO
  go to 399

320 continue
  !! # periodic:
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        DO m=1,meqn
           q(m,i,1-jbc) = q(m,i,my+1-jbc)
        ENDDO
     ENDDO
  ENDDO
  go to 399

330 continue
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

399 continue
  !!
  !!-------------------------------------------------------
  !!     # top boundary:
  !!-------------------------------------------------------
  go to (400,410,420,430) mthbc(4)+1

400 continue
  !! # user-specified boundary conditions go here in place of error output
  write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
  stop
  go to 499

410 continue
  !! # zero-order extrapolation:
  DO jbc=1,mbc
     DO i = 1-mbc, mx+mbc
        DO m=1,meqn
           q(m,i,my+jbc) = q(m,i,my)
        ENDDO
     ENDDO
  ENDDO
  go to 499

420 continue
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

499 continue

  RETURN
END SUBROUTINE clawapck5_bc2
