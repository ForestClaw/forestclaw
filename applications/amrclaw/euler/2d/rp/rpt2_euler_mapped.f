c =========================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,
     &      mx,ql,qr,aux1,aux2, aux3, ilr,
     &      asdq, bmasdq, bpasdq)
c =========================================================
c
c     # solve Riemann problems for the 1D Euler equations using Roe's
c     # approximate Riemann solver.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     # On output, wave contains the waves,
c     #      the speeds,
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp is called with ql = qr = q.
c
c
      implicit none

      integer maxm, meqn, mwaves, mbc, mx, ilr, ixy
      double precision    ql(1-mbc:maxm+mbc, meqn)
      double precision    qr(1-mbc:maxm+mbc, meqn)
      double precision     s(1-mbc:maxm+mbc, mwaves)
      double precision  wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision  asdq(1-mbc:maxm+mbc, meqn)
      double precision  bmasdq(1-mbc:maxm+mbc, meqn)
      double precision  bpasdq(1-mbc:maxm+mbc, meqn)
      double precision  aux1(1-mbc:maxm+mbc,*)
      double precision  aux2(1-mbc:maxm+mbc,*)
      double precision  aux3(1-mbc:maxm+mbc,*)

c     # For Roe solver
      double precision rhol, ul, vl, el, cl, pl
      double precision rhor, ur, vr, er, cr, pr
      double precision u, v, enth, rho, p, e
      double precision rhsqrtl, rhsqrtr, rhsq2
      double precision uvl, uvr, uv(2)


c     # For mappings
      integer meqn2,mwaves2
      parameter(meqn2 = 4, mwaves2 = 3)
      double precision ql_state(meqn2), qr_state(meqn2)
      double precision q_state(meqn2)
      double precision wave_local(meqn2,mwaves2)
      double precision s_local(mwaves2), delta(meqn2)
      double precision speeds(mwaves2,2)
      double precision deltam(meqn2), deltap(meqn2)
      double precision area

c     # for mapping
      double precision rotm(4), rotp(4), uvm_rot(2),uvp_rot(2)
      integer locrot,mcapa,locarea

c     # Miscellaneous
      integer i, j, m, mw, i1, ixy1

c     # Problem parameters
      double precision gamma, gamma1

      double precision dtcom, dxcom, dycom, tcom
      integer icom, jcom

      common /param/  gamma,gamma1
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      call get_aux_locations_t(ixy, mcapa, locrot,locarea)

      do i = 2-mbc,mx+mbc
         i1 = i + ilr - 2

c          do m = 1,meqn
c             ql_state(m) = qr(i-1,m)
c             qr_state(m) = ql(i,m)
c          enddo
c
c          rhol = ql_state(1)
c          rhor = qr_state(1)
c
c          ul = ql_state(2)/rhol
c          ur = qr_state(2)/rhor
c
c          vl = ql_state(3)/rhol
c          vr = qr_state(3)/rhor
c
c          el = ql_state(4)
c          er = qr_state(4)
c
c          uvl = ul*ul + vl*vl
c          uvr = ur*ur + vr*vr
c          pl = gamma1*(el - 0.5d0*rhol*uvl)
c          pr = gamma1*(er - 0.5d0*rhor*uvr)
c
c c        # Get Roe averaged values
c          rhsqrtl = sqrt(rhol)
c          rhsqrtr = sqrt(rhor)
c          rhsq2 = rhsqrtl + rhsqrtr


c         uv(1) = (ul*rhsqrtl + ur*rhsqrtr) / rhsq2
c         uv(2) = (vl*rhsqrtl + vr*rhsqrtr) / rhsq2
c         enth = ((el + pl)/rhsqrtl + (er + pr)/rhsqrtr)/rhsq2

         do m = 1,meqn
            if (ilr .eq. 1) then
               q_state(m) = ql(i,m)
            else
               q_state(m) = qr(i-1,m)
            endif
         enddo

         rho = q_state(1)
         uv(1) = q_state(2)/rho
         uv(2) = q_state(3)/rho
         e = q_state(4)
         p = gamma1*(e - 0.5d0*rho*(uv(1)**2 + uv(2)**2))
         enth = (e + p)/rho

         do j = 1,2
            uvm_rot(j) = uv(j)
            uvp_rot(j) = uv(j)
            rotm(j) = aux2(i1,locrot+j-1)
            rotp(j) = aux3(i1,locrot+j-1)
         enddo
         call compute_tangent(rotm)
         call compute_tangent(rotp)

         call rotate2(rotm,uvm_rot)
         call rotate2(rotp,uvp_rot)

c        # Now we have to fix up asdq(i,5) :
         do m = 1,meqn
            deltap(m) = asdq(i,m)
            deltam(m) = asdq(i,m)
         enddo
         call rotate2(rotm,deltam(2))
         call rotate2(rotp,deltap(2))

c        # ------------------------------------------
c        # Solve for minus side
c        # ------------------------------------------
         ixy1 = 1
         call roe_solver(ixy1,uvm_rot,enth,deltam,
     &         wave_local,s_local)

         area = aux2(i1,locarea)
         do mw = 1,mwaves
            call rotate2_tr(rotm,wave_local(2,mw))
            speeds(mw,1) = area*min(s_local(mw),0.d0)
         enddo

         do m = 1,meqn
            bmasdq(i,m) = 0.d0
            do mw = 1,mwaves
               bmasdq(i,m) = bmasdq(i,m) + speeds(mw,1)*wave_local(m,mw)
            enddo
         enddo

c        # ------------------------------------------
c        # Solve for plus side
c        # ------------------------------------------
         ixy1 = 1
         call roe_solver(ixy1,uvp_rot,enth,deltap,
     &         wave_local,s_local)

         area = aux3(i1,locarea)
         do mw = 1,mwaves
            call rotate2_tr(rotp,wave_local(2,mw))
            speeds(mw,2) = area*max(s_local(mw),0.d0)
         enddo

         do m = 1,meqn
            bpasdq(i,m) = 0.d0
            do mw = 1,mwaves
               bpasdq(i,m) = bpasdq(i,m) + speeds(mw,2)*wave_local(m,mw)
            enddo
         enddo

      enddo !! end of i loop

      return
      end
