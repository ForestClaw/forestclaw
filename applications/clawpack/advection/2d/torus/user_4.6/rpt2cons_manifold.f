      subroutine rpt2cons_manifold(ixy,maxm,meqn,mwaves,mbc,mx,
     &                       ql,qr,aux1,aux2,aux3,imp,asdq,
     &                       bmasdq,bpasdq)
      
      implicit none

      integer ixy, maxm, meqn,mwaves,mbc,mx,maux,imp

      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision   aux1(1-mbc:maxm+mbc, *)
      double precision   aux2(1-mbc:maxm+mbc, *)
      double precision   aux3(1-mbc:maxm+mbc, *)


      integer i, i1, k, idir
      double precision vvecl(3), vvecr(3), vhat, nv(3)
      double precision vrrot, vlrot, g, vhat1

c     # ixy = 1 --> idir = 1
c     # ixy = 2 --> idir = 0
      idir = 2-ixy

      do i = 2-mbc, mx+mbc
          i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq

c         # Lower faces - cell centered velocities
          vvecl(1) = aux1(i1,4)    
          vvecl(2) = aux1(i1,5)
          vvecl(3) = aux1(i1,6)

          vvecr(1) = aux2(i1,4)    
          vvecr(2) = aux2(i1,5)
          vvecr(3) = aux2(i1,6)

c         # x-edge lengths (7)
c         # y-edge lengths (8)
          g = aux2(i1,7+idir)

c         # Get scaled edge normals         
          do k = 1,3
              if (idir .eq. 0) then
                  nv(k) = g*aux2(i1,9+k-1)
              else
                  nv(k) = g*aux2(i1,12+k-1)
              endif
          enddo

          vrrot = vvecr(1)*nv(1) + vvecr(2)*nv(2) + vvecr(3)*nv(3)
          vlrot = vvecl(1)*nv(1) + vvecl(2)*nv(2) + vvecl(3)*nv(3)

          vhat = (vrrot + vlrot)/2.0

          bmasdq(i,1) = min(vhat,0.d0)*asdq(i,1)


c         # Upper faces - cell centered velocities
          vvecl(1) = aux2(i1,4)    
          vvecl(2) = aux2(i1,5)
          vvecl(3) = aux2(i1,6)

          vvecr(1) = aux3(i1,4)    
          vvecr(2) = aux3(i1,5)
          vvecr(3) = aux3(i1,6)

c         # x-edge lengths (7)
c         # y-edge lengths (8)
          g = aux3(i1,7+idir)

c         # Get scaled edge normals         
          do k = 1,3
              if (idir .eq. 0) then
                  nv(k) = g*aux3(i1,9+k-1)
              else
                  nv(k) = g*aux3(i1,12+k-1)
              endif
          enddo

          vrrot = vvecr(1)*nv(1) + vvecr(2)*nv(2) + vvecr(3)*nv(3)
          vlrot = vvecl(1)*nv(1) + vvecl(2)*nv(2) + vvecl(3)*nv(3)

          vhat = (vrrot + vlrot)/2.0

          bpasdq(i,1) = max(vhat,0.d0)*asdq(i,1)

      enddo


      return
      end
