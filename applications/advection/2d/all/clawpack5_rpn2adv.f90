SUBROUTINE clawpack5_rpn2adv(ixy,maxm,meqn,mwaves,maux,mbc, &
     mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
  IMPLICIT NONE

  INTEGER ixy,maxm,meqn,mwaves,maux,mbc,mx
  DOUBLE PRECISION wave(meqn, mwaves, 1-mbc:maxm+mbc)
  DOUBLE PRECISION    s(mwaves, 1-mbc:maxm+mbc)
  DOUBLE PRECISION   ql(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION   qr(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION apdq(meqn,1-mbc:maxm+mbc)
  DOUBLE PRECISION amdq(meqn,1-mbc:maxm+mbc)
  DOUBLE PRECISION auxl(maux,1-mbc:maxm+mbc)
  DOUBLE PRECISION auxr(maux,1-mbc:maxm+mbc)

  INTEGER iface, i, i1

  IF (meqn .NE. 1) THEN
     WRITE(6,*) 'clawpack5_rpn2 : meqn should be equal to 1'
     STOP
  ENDIF

  iface = ixy
  DO i = 2-mbc, mx+mbc
     wave(1,1,i) = ql(1,i) - qr(1,i-1)
     s(1,i) = auxl(ixy,i)
     amdq(1,i) = dmin1(auxl(iface,i), 0.d0) * wave(1,1,i)
     apdq(1,i) = dmax1(auxl(iface,i), 0.d0) * wave(1,1,i)
  END DO

  RETURN
END SUBROUTINE clawpack5_rpn2adv
