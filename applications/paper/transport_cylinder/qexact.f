c     # ---------------------------------------------------
c     # QEXACT - returns exact solution to transport problem
c     #          
c     #                q_t + \div (u(x,y,t) q) = 0
c     # 
c     # Format : 
c     #                qe = QEXACT(x,y,t,flow_flag)
c     #
c     # Input arguments : 
c     #
c     # (x,y,t)     : Specifies the time and location where the exact
c     #               solution is desired.
c     # flow_flag   : Set to 0 for divergent-free flow fields and set to
c     #               1 for divergent flow fields.
c     #
c     # 
c     # This routine requires three auxiliary functions : 
c     # 
c     #    call  velocity_components(x,y,t,u)
c     #
c     # which should compute velocity components (u(1),u(2)) in 
c     # terms of spherical basis functions (not normalized)
c     #
c     # The user must also supply a divergence function
c     #
c     #         divu = map_divergence(x,y,t)  
c     #
c     # that computes the divergence of the velocity field. 
c     #
c     # Finally, the user needs to supply a function that computes the
c     # initial conditions
c     #
c     #              q0 = q0_init(x,y) 
c     # 
c     # The solution proceeds in two steps. If we are solving an 
c     # divergent-free problem, we only need to evolve (x(t), y(t)) 
c     # backwards to the initial time.
c     #
c     # If we have a compressible problem, we first evolve (x(t),y(t))
c     # back to a starting location (x0,y0), and  then evolve x,y, and q
c     # from  the starting location back to final position.
c     # back to final position.
c     # 
c     # --------------------------------------------------------

      double precision function qexact(x,y,tfinal,flow_flag)
      implicit none

      external qexact_rhs
      external solout

      double precision x,y,tfinal
      integer flow_flag

      double precision xc0, yc0, t0
      double precision q0_init, q0

      double precision sigma(4), rtol, atol
      integer itol, iout

      integer Nmax, lwork,nrdens, liwork
      parameter(Nmax=4, nrdens=0)
      parameter(lwork=8*Nmax+5*nrdens+21)
      parameter(liwork=nrdens+21)

      double precision work(lwork), rpar(1)
      integer iwork(liwork), ipar(2), idid

      double precision tol

      logical evolve_q

      integer i

c     # ------------------------------------------
c     # Numerical parameters
c     # ------------------------------------------
      itol = 0
      rtol = 1.d-14
      atol = 1.d-14
      iout = 0

      do i = 1,20
          work(i) = 0
          iwork(i) = 0
      enddo

c     # Evolve from t=t0 to t=tfinal
      t0 = 0

c     # Initial conditions for ODE
      sigma(1) = x
      sigma(2) = y

c     # Needed for backward trace      
      rpar(1) = tfinal

c     # Tracing backwards      
      ipar(1) = 0

c     # This traces the velocity field back to the origin.
      call dopri5(2,qexact_rhs,t0,sigma,tfinal,
     &            rtol,atol,itol,
     &            solout,iout, work,lwork,iwork,liwork,
     &            rpar,ipar,idid)

      if (idid .ne. 1) then
          write(6,*) 'DOPRI5 : idid .ne. 1'
          stop
      endif

c     # Initial position in [0,1]x[0,1]
      xc0 = sigma(1)
      yc0 = sigma(2)

c     # Get initial condition in computational coordinates (xc0,yc0)
      q0 =  q0_init(xc0,yc0)

c     # Evolve q along characteristics for divergent case     
      evolve_q = flow_flag == 1
      if (evolve_q) then
c         # Variable coefficient case        
c         # We now need to evolve q along with (x,y), starting from
c         # from (xc0,yc0)
          sigma(1) = xc0
          sigma(2) = yc0
          sigma(3) = q0

          do i = 1,20
              work(i) = 0
              iwork(i) = 0
          enddo

c         # Tracing forwards
          ipar(1) = 1

          t0 = 0
          call dopri5(3,qexact_rhs,t0,sigma,tfinal,
     &                rtol,atol,itol,
     &                solout,iout, work,lwork,iwork,liwork,
     &                rpar,ipar,idid)

          if (idid .ne. 1) then
              write(6,*) 'DOPRI5 : idid .ne. 1'
              stop
          endif

          tol = atol
          if (abs(sigma(1)-x) .gt. tol) then
              write(6,*) 'qexact.f : Did not evolve x correctly'
              write(6,100) xc0, x, sigma(1), abs(x-sigma(1))
              write(6,*)' '
              stop
          endif
          if (abs(sigma(2)-y) .gt. tol) then
              write(6,*) 'qexact.f : Did not evolve y correctly'
              write(6,100) xc0,x, sigma(2), abs(y-sigma(2))
              write(6,*)' '
              stop
          endif
100       format(3F24.16, E24.4)

          qexact = sigma(3)
      else
          qexact = q0
      endif

      end

c     # ---------------------------------------------------------------      
      subroutine solout(nr,xold,x,y,n,con,icomp,nd,rpar,ipar,irtrn)
      dimension y(n),con(5*nd),icomp(nd)

c     # Dummy routine

      end


      subroutine qexact_rhs(n,t,sigma,f,rpar,ipar)
      implicit none

      integer n, ipar(1)
      double precision t, sigma(n), f(n), rpar(1)
      double precision x,y, u(2), tt, tfinal
      double precision q, divu, map_divergence
      integer idir

      idir = ipar(1)

      x = sigma(1)
      y = sigma(2)

      if (idir .eq. 0) then
          tfinal = rpar(1)
          tt = tfinal - t
      else
          tt = t
      endif

      call velocity_components(x,y,tt,u)

      if (idir .eq. 0) then
c         # We are tracing back to initial position
          f(1) = -u(1)
          f(2) = -u(2)
      else
c         # We are evolving q along characteristics starting 
c         # from initial position found above
          q = sigma(3)
          divu = map_divergence(x,y,t)

          f(1) = u(1)
          f(2) = u(2)
          f(3) = -divu*q   
      endif
      end
