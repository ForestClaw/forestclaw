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

  INTEGER iface, i, i1, mq

  iface = 3-ixy  !#  = 1 if ixy=2  or  = 2 if ixy=1
  DO i = 2-mbc,mx+mbc
     i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
     do mq = 1,meqn
        bmasdq(mq,i) = min(aux2(iface,i1), 0.d0) * asdq(mq,i)
        bpasdq(mq,i) = max(aux3(iface,i1), 0.d0) * asdq(mq,i)
     end do
  END DO

  RETURN
END SUBROUTINE clawpack5_rpt2adv
