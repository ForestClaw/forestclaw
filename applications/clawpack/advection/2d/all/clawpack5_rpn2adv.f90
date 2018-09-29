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

    INTEGER iface, i, i1, mq

    iface = ixy
    DO i = 2-mbc, mx+mbc
        do mq = 1,meqn
            wave(mq,1,i) = ql(mq,i) - qr(mq,i-1)
        enddo
        s(1,i) = auxl(ixy,i)
        DO mq = 1,meqn
            amdq(mq,i) = dmin1(auxl(iface,i), 0.d0) * wave(mq,1,i)
            apdq(mq,i) = dmax1(auxl(iface,i), 0.d0) * wave(mq,1,i)
        ENDDO
    END DO

    RETURN
END SUBROUTINE clawpack5_rpn2adv
