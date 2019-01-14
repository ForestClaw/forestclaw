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
      double precision vrrot, vlrot, g

c     # ixy = 1 --> idir = 1
c     # ixy = 2 --> idir = 0
      idir = 2-ixy

      do i = 2-mbc, mx+mbc
          i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq

c         # -----------------------------------------
c         # Lower faces - cell centered velocities
c         # -----------------------------------------
c          vvecl(1) = aux1(i1,4)    
c          vvecl(2) = aux1(i1,5)
c          vvecl(3) = aux1(i1,6)

c         vvecr(1) = aux2(i1,4)    
c         vvecr(2) = aux2(i1,5)
c         vvecr(3) = aux2(i1,6)

c         # x-edge lengths (7)
c         # y-edge lengths (8)
c          g = aux2(i1,7+idir)

c         # Get scaled edge normals         
c          do k = 1,3
c              if (idir .eq. 0) then
c                  nv(k) = g*aux2(i1,9+k-1)
c              else
c                  nv(k) = g*aux2(i1,12+k-1)
c              endif
c          enddo

c          vrrot = vvecr(1)*nv(1) + vvecr(2)*nv(2) + vvecr(3)*nv(3)
c          vlrot = vvecl(1)*nv(1) + vvecl(2)*nv(2) + vvecl(3)*nv(3)

           
          g = aux2(i1,6+idir)

          vrrot = aux2(i,  2 + 2*idir)   !! Left edge of right cell
          vlrot = aux2(i-1,3 + 2*idir)   !! Right edge of left cell

          vhat = g*(vrrot + vlrot)/2.0

          bmasdq(i,1) = min(vhat,0.d0)*asdq(i,1)

c         # -----------------------------------------
c         # Upper faces - cell centered velocities
c         # -----------------------------------------
c          vvecl(1) = aux2(i1,4)    
c          vvecl(2) = aux2(i1,5)
c          vvecl(3) = aux2(i1,6)

c          vvecr(1) = aux3(i1,4)    
c          vvecr(2) = aux3(i1,5)
c          vvecr(3) = aux3(i1,6)

c         # x-edge lengths (7)
c         # y-edge lengths (8)
c          g = aux3(i1,7+idir)

c         # Get scaled edge normals         
c          do k = 1,3
c              if (idir .eq. 0) then
c                  nv(k) = g*aux3(i1,9+k-1)
c              else
c                  nv(k) = g*aux3(i1,12+k-1)
c              endif
c          enddo

c          vrrot = vvecr(1)*nv(1) + vvecr(2)*nv(2) + vvecr(3)*nv(3)
c          vlrot = vvecl(1)*nv(1) + vvecl(2)*nv(2) + vvecl(3)*nv(3)

          g = aux3(i1,6+idir)

          vrrot = aux3(i,  2 + 2*idir)   !! Left edge of right cell
          vlrot = aux3(i-1,3 + 2*idir)   !! Right edge of left cell

          vhat = g*(vrrot + vlrot)/2.0

          bpasdq(i,1) = max(vhat,0.d0)*asdq(i,1)

      enddo


      return
      end
