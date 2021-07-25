c     # ------------------------------------------------------------
c     # Streamfunction, velocity components and derivatives
c     # 
c     # ------------------------------------------------------------
      double precision function torus_psi(x,y)
      implicit none

      double precision x, y

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      double precision revs_per_s
      common /stream_comm/ revs_per_s

      double precision psi

      if (beta .ne. 0) then
          write(6,'(A,A)') 'psi (psi.f) : Streamfunction only works ',
     &         'for beta == 0'
          stop 'error'
      endif

      psi = pi2*revs_per_s*alpha*(pi2*y + alpha*sin(pi2*y))

      torus_psi = psi

      end

c     # ------------------------------------------------------------
c     # Edge : u = \nabla cross \Psi - curl \Psi
c     # NOTE : This depends on computational coordinates, not 
c     # physical coordinates.
c     # ------------------------------------------------------------
      subroutine torus_edge_velocity(x1,y1,x2,y2,ds,vn)
      implicit none

      double precision x1,y1,x2,y2, ds, vn, torus_psi

      vn = (torus_psi(x1,y1) - torus_psi(x2,y2))/ds

      end

