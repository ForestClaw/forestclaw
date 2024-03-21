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

    INTEGER iface, i, mq, mw


    double precision delta(meqn)

    iface = ixy
    DO i = 2-mbc, mx+mbc
        do mq = 1,meqn
            do mw = 1,mwaves
                wave(mq,mw,i) = 0
            end do
        end do

        do mq = 1,meqn
            delta(mq) = ql(mq,i) - qr(mq,i-1)
            do mw = 1,mwaves
                wave(mq,mw,i) = delta(mq)
            end do
        end do

        do mw = 1,mwaves
            s(mw,i) = auxl(iface,i)
        end do

        do mq = 1,meqn
            amdq(mq,i) = 0
            apdq(mq,i) = 0
            do mw = 1,mwaves
               amdq(mq,i) = amdq(mq,i) + min(s(mw,i), 0.d0) * wave(mq,mw,i)
               apdq(mq,i) = apdq(mq,i) + max(s(mw,i), 0.d0) * wave(mq,mw,i)
            end do
        end do
    END DO

    RETURN
END SUBROUTINE clawpack5_rpn2adv
