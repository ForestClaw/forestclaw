      subroutine rpn2(ixy,maxm, meqn,mwaves,mbc,mx, ql,qr,
     &      auxl,auxr,wave,s,amdq,apdq)
      implicit none

      integer ixy, maxm, meqn, mwaves, mbc, mx, my
      double precision    ql(1-mbc:maxm+mbc, meqn)
      double precision    qr(1-mbc:maxm+mbc, meqn)
      double precision     s(1-mbc:maxm+mbc, mwaves)
      double precision  wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision  amdq(1-mbc:maxm+mbc, meqn)
      double precision  apdq(1-mbc:maxm+mbc, meqn)
      double precision  auxl(1-mbc:maxm+mbc,*)
      double precision  auxr(1-mbc:maxm+mbc,*)

      double precision gamma, gamma1
      double precision ql_state(4), qr_state(4)
      double precision rhol, ul, el, cl, pl, vl
      double precision rhor, ur, er, cr, pr, vr
      double precision rhsqrtl, rhsqrtr, rhsq2
      double precision u, enth, delta(4), rho, v
      double precision wave_local(4,3), s_local(3),uv(2)
      double precision speeds(3,2), u2v2l, u2v2r
      integer i, m, mw, mu, mv, ixy1
      logical efix

      integer mcapa,locrot, locarea
      double precision rot(4), area

      common /param/  gamma,gamma1

      data efix /.true./


      call get_aux_locations_n(ixy,mcapa,locrot,locarea)

      do i = 2-mbc,mx+mbc

         rot(1) = auxl(i,locrot)
         rot(2) = auxl(i,locrot+1)
         call compute_tangent(rot)

         do m = 1,meqn
            ql_state(m) = qr(i-1,m)
            qr_state(m) = ql(i,m)
         enddo
         call rotate2(rot,ql_state(2))
         call rotate2(rot,qr_state(2))

         rhol = ql_state(1)
         rhor = qr_state(1)

         ul = ql_state(2)/rhol
         ur = qr_state(2)/rhor

         vl = ql_state(3)/rhol
         vr = qr_state(3)/rhor

         el = ql_state(4)
         er = qr_state(4)

         u2v2l = ul*ul + vl*vl
         u2v2r = ur*ur + vr*vr
         pl = gamma1*(el - 0.5d0*rhol*u2v2l)
         pr = gamma1*(er - 0.5d0*rhor*u2v2r)

c        # Get Roe averaged values
         rhsqrtl = sqrt(rhol)
         rhsqrtr = sqrt(rhor)
         rhsq2 = rhsqrtl + rhsqrtr

         uv(1) = (ul*rhsqrtl + ur*rhsqrtr) / rhsq2
         uv(2) = (vl*rhsqrtl + vr*rhsqrtr) / rhsq2
         enth = ((el + pl)/rhsqrtl + (er + pr)/rhsqrtr)/rhsq2

         do m = 1,meqn
            delta(m) = qr_state(m) - ql_state(m)
         enddo

         ixy1 = 1
         call roe_solver(ixy1,uv,enth,delta,wave_local,s_local)

         do mw = 1,mwaves
            speeds(mw,1) = min(s_local(mw),0.d0)
            speeds(mw,2) = max(s_local(mw),0.d0)
         enddo

         if (efix) then
c           # This modifies the speeds, but we will still have
c           # s(mw) = speeds(mw,1) + speeds(mw,2)
            cl = sqrt(gamma*pl/rhol)
            cr = sqrt(gamma*pr/rhor)
            call apply_entropy_fix(ql_state,qr_state,cl,cr,
     &            wave_local, speeds)
         endif

         area = auxl(i,locarea)
         do mw = 1,mwaves
            call rotate2_tr(rot,wave_local(2,mw))
            speeds(mw,1) = area*speeds(mw,1)
            speeds(mw,2) = area*speeds(mw,2)
         enddo

         do m = 1,meqn
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do mw = 1,mwaves
               wave(i,m,mw) = wave_local(m,mw)
               s(i,mw) = speeds(mw,1) + speeds(mw,2)
               amdq(i,m) = amdq(i,m) + speeds(mw,1)*wave(i,m,mw)
               apdq(i,m) = apdq(i,m) + speeds(mw,2)*wave(i,m,mw)
            enddo
         enddo
      enddo


      return
      end

      subroutine apply_entropy_fix(ql,qr,cl, cr,wave_local,speeds)
      implicit none

      double precision ql(4), qr(4),wave_local(4,3)
      double precision speeds(4,2), s1, s2, s3, s4
      double precision sl, sml, smr, sr, qml(4),qmr(4)
      double precision ul, cl, ur, cr, pml, pmr, cml,cmr
      double precision sfract
      double precision gamma, gamma1
      logical trans1

      integer m

      common /param/ gamma, gamma1

      s1 = speeds(1,1) + speeds(1,2)
      s2 = speeds(2,1) + speeds(2,2)
      s3 = speeds(3,1) + speeds(3,2)

      do m = 1,4
         qml(m) = ql(m) + wave_local(m,1)
      enddo
      sl = ql(2)/ql(1) - cl
      pml = gamma1*(qml(4) - 0.5d0*(qml(2)**2/qml(1) +
     &      qml(3)**2/qml(1)))
      sml = qml(2)/qml(1) - sqrt(gamma*pml/qml(1))

c     # Check the 1-wave
      trans1 = .false.
      if (sl .lt. 0 .and. sml .gt. 0) then
c        # apply transonic entropy fix
         trans1 = .true.
         s1 = (sl + sml)/2.d0
         sfract = (sml - s1)/(sml - sl)
         speeds(1,1) = sfract*sl
         speeds(1,2) = (1-sfract)*sml
      endif

c     # If the 1-wave is transonic,then we are done...
c     # Otherwise, we have to check the 3-wave
      if (.not. trans1) then
         do m = 1,4
            qmr(m) = qr(m) - wave_local(m,3)
         enddo
         sr = qr(2)/qr(1) + cr
         pmr = gamma1*(qmr(4) - 0.5d0*(qmr(2)**2/qmr(1) +
     &         qmr(3)**2/qmr(1)))
         smr = qmr(2)/qmr(1) + sqrt(gamma*pmr/qmr(1))
         if (smr .lt. 0 .and. sr .gt. 0) then
c           # apply transonic entropy fix
            s3 = (smr + sr)/2.d0
            sfract = (sr - s3)/(sr - smr)
            speeds(3,1) = sfract*smr
            speeds(3,2) = (1-sfract)*sr
         endif
      endif

      end
