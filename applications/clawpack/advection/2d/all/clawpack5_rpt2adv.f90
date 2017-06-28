SUBROUTINE clawpack5_rpt2adv(ixy,imp,maxm,meqn,mwaves,maux,&
     mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
  IMPLICIT NONE

  INTEGER ixy,imp,maxm,meqn,mwaves,maux,mbc,mx

  DOUBLE PRECISION     ql(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION     qr(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION   asdq(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION bmasdq(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION bpasdq(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION   aux1(maux, 1-mbc:maxm+mbc)
  DOUBLE PRECISION   aux2(maux, 1-mbc:maxm+mbc)
  DOUBLE PRECISION   aux3(maux, 1-mbc:maxm+mbc)

  INTEGER iface, i, i1

  IF (meqn .NE. 1) THEN
     WRITE(6,*) 'clawpack5_rpn2 : meqn should be equal to 1'
     STOP
  ENDIF


  iface = 3-ixy  !#  = 1 if ixy=2  or  = 2 if ixy=1
  DO i = 2-mbc,mx+mbc
     i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
     bmasdq(1,i) = dmin1(aux2(iface,i1), 0.d0) * asdq(1,i)
     bpasdq(1,i) = dmax1(aux3(iface,i1), 0.d0) * asdq(1,i)
  END DO

  RETURN
END SUBROUTINE clawpack5_rpt2adv
