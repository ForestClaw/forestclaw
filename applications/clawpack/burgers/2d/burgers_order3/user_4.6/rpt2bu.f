c
cc     =====================================================
c      subroutine clawpack46_rpt2adv(ixy,maxm,meqn,mwaves,mbc,mx,
c     &                  ql,qr,aux1,aux2,aux3, ilr,asdq,bmasdq,bpasdq)
cc     =====================================================

c
c     =====================================================
      subroutine clawpack46_rpt2bu(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  q,ql,aux1,aux2,aux3,
     &                  imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for 2D Burgers' equation
c
c     # Split asdq into eigenvectors of Roe matrix B.
c     # For the scalar equation, this simply amounts to computing the
c     # transverse wave speed from the opposite Riemann problem.
c
      dimension     q(1-mbc:maxm+mbc, meqn)
      dimension   asdq(1-mbc:maxm+mbc,meqn)
      dimension bmasdq(1-mbc:maxm+mbc, meqn)
      dimension bpasdq(1-mbc:maxm+mbc, meqn)



c
c     # x- and y- Riemann problems are identical, so it doesn't matter if
c     # ixy=1 or 2.
c
          do 10 i = 4-mbc, mx+mbc-3
            
c            sb = (7d0*(q(i,1)+q(i-1,1))-q(i-2,1)-q(i+1,1))/12d0
             sb = 0.5*(q(i,1)+q(i-1,1))
             bmasdq(i,1) = dmin1(sb, 0.d0) * asdq(i,1)
             bpasdq(i,1) = dmax1(sb, 0.d0) * asdq(i,1)
   10        continue

      return
      end
