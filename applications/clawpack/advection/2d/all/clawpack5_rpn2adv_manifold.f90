SUBROUTINE clawpack5_rpn2adv_manifold(ixy,maxm,meqn, &
     mwaves,maux,mbc, mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
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

  INTEGER iface, i, i1, m

  IF (mwaves .NE. 1) THEN
     WRITE(6,*) 'clawpack5_rpn2 : mwaves should be equal to 1'
     STOP
  ENDIF

  iface = ixy
  DO i = 2-mbc, mx+mbc
     DO m = 1,meqn
        wave(m,1,i) = ql(m,i) - qr(m,i-1)
     ENDDO

     !! # Only one wave
     s(1,i) = auxl(iface+1,i)

     DO m = 1,meqn
        amdq(m,i) = MIN(s(1,i), 0.d0) * wave(m,1,i)
        apdq(m,i) = MAX(s(1,i), 0.d0) * wave(m,1,i)
     ENDDO
  END DO

  RETURN
END SUBROUTINE clawpack5_rpn2adv_manifold
