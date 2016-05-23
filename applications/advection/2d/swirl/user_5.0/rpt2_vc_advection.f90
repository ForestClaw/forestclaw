! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
    implicit double precision(a-h,o-z)

!     # Riemann solver in the transverse direction for the advection equation.

    dimension     ql(meqn, 1-mbc:maxm+mbc)
    dimension     qr(meqn, 1-mbc:maxm+mbc)
    dimension   asdq(meqn, 1-mbc:maxm+mbc)
    dimension bmasdq(meqn, 1-mbc:maxm+mbc)
    dimension bpasdq(meqn, 1-mbc:maxm+mbc)
    dimension   aux1(maux, 1-mbc:maxm+mbc)
    dimension   aux2(maux, 1-mbc:maxm+mbc)
    dimension   aux3(maux, 1-mbc:maxm+mbc)

    kv = 3-ixy  !#  = 1 if ixy=2  or  = 2 if ixy=1
    do 10 i=2-mbc,mx+mbc
        i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
        bmasdq(1,i) = dmin1(aux2(kv,i1), 0.d0) * asdq(1,i)
        bpasdq(1,i) = dmax1(aux3(kv,i1), 0.d0) * asdq(1,i)
    10 END DO

    return
    end subroutine rpt2
