      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &      ql,qr,aux1,aux2,aux3,
     &      imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for the shallow water
c     # equations  on a quadrilateral grid.
c
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # Uses Roe averages and other quantities which were
c     # computed in rpn2sh and stored in the common block comroe.
c
      dimension     ql(1-mbc:maxm+mbc, meqn)
      dimension     qr(1-mbc:maxm+mbc, meqn)
      dimension   asdq(1-mbc:maxm+mbc, meqn)
      dimension bmasdq(1-mbc:maxm+mbc, meqn)
      dimension bpasdq(1-mbc:maxm+mbc, meqn)
      dimension   aux1(1-mbc:maxm+mbc, *)
      dimension   aux2(1-mbc:maxm+mbc, *)
      dimension   aux3(1-mbc:maxm+mbc, *)
c
      parameter (maxm2 = 1002)          !# assumes at most 1000x1000 grid with mbc=2
      dimension u(-1:maxm2),v(-1:maxm2),a(-1:maxm2),h(-1:maxm2)
      dimension wave(-1:maxm2, 4, 3)
      dimension    s(-1:maxm2, 3)
      dimension enx(-1:maxm2), eny(-1:maxm2), enz(-1:maxm2)
      dimension etx(-1:maxm2), ety(-1:maxm2), etz(-1:maxm2)
      dimension gamma(-1:maxm2)

      dimension delta(4)
      common /sw/  g
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c

      if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
         write(6,*) 'need to increase maxm2 in rpB'
         stop
      endif

      if(ixy.eq.1) then
         dx = dxcom
      else
         dx = dycom
      endif
c
c
      if(ixy.eq.1) then
         ioff = 7
      else
         ioff = 1
      endif

c
c     # imp is used to flag whether wave is going to left or right,
c     # since states and grid are different on each side
c
      if (imp.eq.1) then
c        # asdq = amdq, moving to left
         ix1 = 2-mbc
         ixm1 = mx+mbc
      else
c        # asdq = apdq, moving to right
         ix1 = 1-mbc
         ixm1 = mx+mbc
      endif
c
c     --------------
c     # up-going:
c     --------------
c

c     # determine rotation matrix for interface above cell, using
c     aux3
c     [ alf  beta ]
c     [-beta  alf ]
c
      do i=ix1,ixm1
c
         if (imp.eq.1) then
            i1 = i-1
         else
            i1 = i
         endif
c
         enx(i) =   aux3(i1,ioff+1)
         eny(i) =   aux3(i1,ioff+2)
         enz(i) =   aux3(i1,ioff+3)
         etx(i) =   aux3(i1,ioff+4)
         ety(i) =   aux3(i1,ioff+5)
         etz(i) =   aux3(i1,ioff+6)
         gamma(i) = dsqrt(etx(i)**2 + ety(i)**2 + etz(i)**2)
         etx(i) =   etx(i) / gamma(i)
         ety(i) =   ety(i) / gamma(i)
         etz(i) =   etz(i) / gamma(i)

         h(i) = ql(i1,1)
         u(i) = (enx(i)*ql(i1,2)+eny(i)*ql(i1,3)+enz(i)*ql(i1,4))
     &         / h(i)
         v(i) = (etx(i)*ql(i1,2)+ety(i)*ql(i1,3)+etz(i)*ql(i1,4))
     &         / h(i)
         a(i) = dsqrt(g*h(i))
      enddo
c
c
c
c     # now split asdq into waves:
c
c     # find a1 thru a3, the coefficients of the 3 eigenvectors:
      do 20 i = ix1,ixm1
         delta(1) = asdq(i,1)
         delta(2) = enx(i)*asdq(i,2)+eny(i)*asdq(i,3)+enz(i)*asdq(i,4)
         delta(3) = etx(i)*asdq(i,2)+ety(i)*asdq(i,3)+etz(i)*asdq(i,4)
         a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
         a2 = -v(i)*delta(1) + delta(3)
         a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))
c
c        # Compute the waves.
c
         wave(i,1,1) = a1
         wave(i,2,1) = enx(i)*a1*(u(i)-a(i)) + etx(i)*a1*v(i)
         wave(i,3,1) = eny(i)*a1*(u(i)-a(i)) + ety(i)*a1*v(i)
         wave(i,4,1) = enz(i)*a1*(u(i)-a(i)) + etz(i)*a1*v(i)
         s(i,1) = (u(i)-a(i))*gamma(i)/dx
c
         wave(i,1,2) = 0.0d0
         wave(i,2,2) = etx(i)*a2
         wave(i,3,2) = ety(i)*a2
         wave(i,4,2) = etz(i)*a2
         s(i,2) = u(i) * gamma(i)/dx
c
         wave(i,1,3) = a3
         wave(i,2,3) = enx(i)*a3*(u(i)+a(i)) + etx(i)*a3*v(i)
         wave(i,3,3) = eny(i)*a3*(u(i)+a(i)) + ety(i)*a3*v(i)
         wave(i,4,3) = enz(i)*a3*(u(i)+a(i)) + etz(i)*a3*v(i)
         s(i,3) = (u(i)+a(i)) * gamma(i)/dx
   20 continue
c
c
c     # compute flux difference bpasdq
c     --------------------------------
c
      do 40 m=1,meqn
         do 40 i=ix1,ixm1
	    bpasdq(i,m) = 0.d0
	    do 30 mw=1,mwaves
	       bpasdq(i,m) = bpasdq(i,m)
     &               + dmax1(s(i,mw),0.d0)*wave(i,m,mw)
   30       continue
   40    continue
c
c        # project momentum component of bpasdq to tangent plane:
         do i=ix1, ixm1
            if (imp.eq.1) then
               i1 = i-1
            else
               i1 = i
            endif
            erx = aux3(i1,14)
            ery = aux3(i1,15)
            erz = aux3(i1,16)
            bn = erx*bpasdq(i,2)+ery*bpasdq(i,3)+erz*bpasdq(i,4)
            bpasdq(i,2) = bpasdq(i,2) - bn*erx
            bpasdq(i,3) = bpasdq(i,3) - bn*ery
            bpasdq(i,4) = bpasdq(i,4) - bn*erz

         enddo

c
c
c        --------------
c        # down-going:
c        --------------
c

c        # determine rotation matrix for interface below cell, using
c        aux2
c        [ alf  beta ]
c        [-beta  alf ]
c
         do i=ix1,ixm1
c
            if (imp.eq.1) then
               i1 = i-1
            else
               i1 = i
            endif
c
            enx(i) =   aux2(i1,ioff+1)
            eny(i) =   aux2(i1,ioff+2)
            enz(i) =   aux2(i1,ioff+3)
            etx(i) =   aux2(i1,ioff+4)
            ety(i) =   aux2(i1,ioff+5)
            etz(i) =   aux2(i1,ioff+6)
            gamma(i) = dsqrt(etx(i)**2 + ety(i)**2 + etz(i)**2)
            etx(i) =   etx(i) / gamma(i)
            ety(i) =   ety(i) / gamma(i)
            etz(i) =   etz(i) / gamma(i)
            u(i) = (enx(i)*ql(i1,2)+eny(i)*ql(i1,3)+enz(i)*ql(i1,4))
     &            / h(i)
            v(i) = (etx(i)*ql(i1,2)+ety(i)*ql(i1,3)+etz(i)*ql(i1,4))
     &            / h(i)
         enddo
c
c
c
c        # now split asdq into waves:
c
c        # find a1 thru a3, the coefficients of the 3 eigenvectors:
         do 80 i = ix1,ixm1
            delta(1) = asdq(i,1)
            delta(2) = enx(i)*asdq(i,2)+eny(i)*asdq(i,3)+enz(i)*asdq(i
     &            ,4)
            delta(3) = etx(i)*asdq(i,2)+ety(i)*asdq(i,3)+etz(i)*asdq(i
     &            ,4)

            a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
            a2 = -v(i)*delta(1) + delta(3)
            a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))
c
c           # Compute the waves.
c
            wave(i,1,1) = a1
            wave(i,2,1) = enx(i)*a1*(u(i)-a(i)) + etx(i)*a1*v(i)
            wave(i,3,1) = eny(i)*a1*(u(i)-a(i)) + ety(i)*a1*v(i)
            wave(i,4,1) = enz(i)*a1*(u(i)-a(i)) + etz(i)*a1*v(i)
            s(i,1) = (u(i)-a(i)) * gamma(i)/dx
c
            wave(i,1,2) = 0.0d0
            wave(i,2,2) = etx(i)*a2
            wave(i,3,2) = ety(i)*a2
            wave(i,4,2) = etz(i)*a2
            s(i,2) = u(i) * gamma(i)/dx
c
            wave(i,1,3) = a3
            wave(i,2,3) = enx(i)*a3*(u(i)+a(i)) + etx(i)*a3*v(i)
            wave(i,3,3) = eny(i)*a3*(u(i)+a(i)) + ety(i)*a3*v(i)
            wave(i,4,3) = enz(i)*a3*(u(i)+a(i)) + etz(i)*a3*v(i)
            s(i,3) = (u(i)+a(i)) * gamma(i)/dx
   80    continue
c
c
c        # compute flux difference bmasdq
c        --------------------------------
c
         do 100 m=1,meqn
            do 100 i=ix1,ixm1
               bmasdq(i,m) = 0.d0
               do 90 mw=1,mwaves
                  bmasdq(i,m) = bmasdq(i,m)
     &                  + dmin1(s(i,mw), 0.d0)*wave(i,m,mw)
   90          continue
  100       continue
c           # project momentum component of bmasdq to tangent plane:
            do i=ix1, ixm1
               if (imp.eq.1) then
                  i1 = i-1
               else
                  i1 = i
               endif
               erx = aux1(i1,14)
               ery = aux1(i1,15)
               erz = aux1(i1,16)
               bn = erx*bmasdq(i,2)+ery*bmasdq(i,3)+erz*bmasdq(i,4)
               bmasdq(i,2) = bmasdq(i,2) - bn*erx
               bmasdq(i,3) = bmasdq(i,3) - bn*ery
               bmasdq(i,4) = bmasdq(i,4) - bn*erz

            enddo
c
c

            return
            end
