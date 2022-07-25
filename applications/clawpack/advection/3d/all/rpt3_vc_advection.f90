subroutine clawpack46_rpt3(ixyz,icoor,imp,maxm,meqn,mwaves,maux,mbc, & 
    mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
    ! Riemann solver in the transverse direction for the
    ! advection equations.
    !
    ! On input,
    !
    !    ql,qr is the data along some one-dimensional slice, as in rpn3
    !         This slice is
    !             in the x-direction if ixyz=1,
    !             in the y-direction if ixyz=2, or
    !             in the z-direction if ixyz=3.
    !    asdq is an array of flux differences (A^*\Dq).
    !         asdq(i,:) is the flux difference propagating away from
    !         the interface between cells i-1 and i.
    !    Note that asdq represents B^*\Dq if ixyz=2 or C^*\Dq if ixyz=3.
    !
    !    ixyz indicates the direction of the original Riemann solve,
    !         called the x-like direction in the table below:
    !
    !               x-like direction   y-like direction   z-like direction
    !      ixyz=1:        x                  y                  z
    !      ixyz=2:        y                  z                  x
    !      ixyz=3:        z                  x                  y
    !
    !    icoor indicates direction in which the transverse solve should
    !         be performed.
    !      icoor=2: split in the y-like direction.
    !      icoor=3: split in the z-like direction.
    !
    !    For example,
    !      ixyz=1, icoor=2 means asdq=A^*\Dq, and should be split in y into
    !                        bmasdq = B^-A^*\Dq,
    !                        bpasdq = B^+A^*\Dq.
    !
    !      ixyz=2, icoor=2 means asdq=B^*\Dq, and should be split in z into
    !                        bmasdq = C^-B^*\Dq,
    !                        bpasdq = C^+B^*\Dq.
    !
    !    The parameter imp is generally needed only if aux
    !    arrays are being used, in order to access the appropriate
    !    variable coefficients:
    !
    !    imp = 1 if asdq = A^- \Dq,  the left-going flux difference
    !          2 if asdq = A^+ \Dq, the right-going flux difference
    !
    !    aux2(:,:,2) is a 1d slice of the aux array along the row
    !                 where the data ql, qr lie.   
    !    aux1(:,:,2) and aux3(:,:,2) are neighboring rows in the 
    !                 y-like direction
    !    aux2(:,:,1) and aux2(:,:,3) are neighboring rows in the 
    !                z-like direction

    implicit none

    integer, intent(in) :: ixyz, icoor, maxm, meqn, mwaves, mbc, mx, maux, imp
    double precision, intent(in) :: ql, qr, asdq, aux1, aux2, aux3
    double precision, intent(out) :: bmasdq, bpasdq

    dimension     ql(meqn,1-mbc:maxm+mbc)
    dimension     qr(meqn,1-mbc:maxm+mbc)
    dimension   asdq(meqn,1-mbc:maxm+mbc)
    dimension bmasdq(meqn,1-mbc:maxm+mbc)
    dimension bpasdq(meqn,1-mbc:maxm+mbc)
    dimension   aux1(maux,1-mbc:maxm+mbc,3)
    dimension   aux2(maux,1-mbc:maxm+mbc,3)
    dimension   aux3(maux,1-mbc:maxm+mbc,3)

    integer manifold
    common /com_manifold/ manifold

    integer :: iuvw, i, i1, mcapa

    ! set iuvw = 1 for u, 2 for v, 3 for w component of velocity
    ! depending on transverse direction:
    mcapa = manifold
    iuvw = ixyz + icoor - 1
    if (iuvw.gt.3) iuvw = iuvw-3

    do i=2-mbc,mx+mbc
        i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
        if (icoor == 2) then  !! transverse dir. is y-like direction
            bmasdq(1,i) = dmin1(aux2(iuvw+mcapa,i1,2), 0.d0)*asdq(1,i)
            bpasdq(1,i) = dmax1(aux3(iuvw+mcapa,i1,2), 0.d0)*asdq(1,i)
        else !! icoor == 3  !! transverse dir. is z-like direction
            !! quanities split into cmasdq and cpasdq
            bmasdq(1,i) = dmin1(aux2(iuvw+mcapa,i1,2),0.d0)*asdq(1,i)
            bpasdq(1,i) = dmax1(aux2(iuvw+mcapa,i1,3),0.d0)*asdq(1,i)
        endif
    end do

    return
end subroutine clawpack46_rpt3

