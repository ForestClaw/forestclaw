! =====================================================
SUBROUTINE clawpack46_rpt2_manifold(ixy,maxm,meqn,mwaves, &
     mbc,mx,ql,qr,aux1,aux2,aux3,imp,asdq,bmasdq,bpasdq)
! =====================================================
! Riemann solver in the transverse direction for the acoustics equations
! on general quadrilateral grids.

! Split asdq (= A^* \Delta q, where * = + or -)
! into down-going flux difference bmasdqb (= B^- A^* \Delta q)
! and up-going flux difference bpasdq (= B^+ A^* \Delta q)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: maxm, meqn, mwaves, mbc, mx, imp
  DOUBLE PRECISION, INTENT(in)      :: ql(1-mbc:maxm+mbc, meqn)
  DOUBLE PRECISION, INTENT(in)      :: qr(1-mbc:maxm+mbc, meqn)
  DOUBLE PRECISION, INTENT(out)   :: asdq(1-mbc:maxm+mbc, meqn)
  DOUBLE PRECISION, INTENT(out) :: bmasdq(1-mbc:maxm+mbc, meqn)
  DOUBLE PRECISION, INTENT(out) :: bpasdq(1-mbc:maxm+mbc, meqn)
  DOUBLE PRECISION, INTENT(in)    :: aux1(1-mbc:maxm+mbc, *)
  DOUBLE PRECISION, INTENT(in)    :: aux2(1-mbc:maxm+mbc, *)
  DOUBLE PRECISION, INTENT(in)    :: aux3(1-mbc:maxm+mbc, *)

  INTEGER :: ixy, ixtran, iytran, ilenrat, i1, i, m
  DOUBLE PRECISION :: cm, cp, zm, zz, zp, xtran, ytran, asdqt, a1, a2
  DOUBLE PRECISION :: asdqp, bpasdqp, bpasdqt, bmasdqp, bmasdqt

  IF (ixy.EQ.1) THEN
     ixtran = 4
     iytran = 5
     ilenrat = 6
  ELSE
     ixtran = 1
     iytran = 2
     ilenrat = 3
  ENDIF

  DO i=2-mbc,mx+mbc
     IF (imp.EQ.1) THEN
        i1 = i-1
     ELSE
        i1 = i
     ENDIF

!       pressure component of asdq:
     asdqp = asdq(i,1)
!
!       sound speed in each row of cells:
     cm = aux1(i1,9)
     cp = aux3(i1,9)
!
!       impedances:
     zm = aux1(i1,8)
     zz = aux2(i1,8)
     zp = aux3(i1,8)

!       up-going:
     xtran = aux3(i,ixtran)
     ytran = aux3(i,iytran)

     asdqt = xtran*asdq(i,2) + ytran*asdq(i,3)
     a2 = (asdqp + zz*asdqt) / (zz+zp)
     bpasdqp = a2*zp
     bpasdqt = a2

     bpasdq(i,1) = cp * bpasdqp
     bpasdq(i,2) = cp * xtran*bpasdqt
     bpasdq(i,3) = cp * ytran*bpasdqt

!       down-going:

     xtran = aux2(i,ixtran)
     ytran = aux2(i,iytran)

     asdqt = xtran*asdq(i,2) + ytran*asdq(i,3)
     a1 = (-asdqp + zz*asdqt) / (zm + zz)
     bmasdqp = -a1*zm
     bmasdqt = a1

     bmasdq(i,1) = -cm * bmasdqp
     bmasdq(i,2) = -cm * xtran*bmasdqt
     bmasdq(i,3) = -cm * ytran*bmasdqt

!       scale by ratio of lengths:
     DO m=1,3
        bmasdq(i,m) = bmasdq(i,m)*aux2(i,ilenrat)
        bpasdq(i,m) = bpasdq(i,m)*aux3(i,ilenrat)
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE clawpack46_rpt2_manifold
