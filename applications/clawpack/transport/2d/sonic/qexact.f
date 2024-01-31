c     # ---------------------------------------------------
c     # Trace back solution to initial position     
c     # Input arguments should be the coordinates
c     # for the canonical swirl mapping in terms of 
c     # (x,y)
c     #       
c     #   T = (R*cos(2*pi*x), R*sin(2*pi*x), r*sin(2*pi*y))
c     #
c     #        r = alpha*(1 + beta*sin(2*pi*x))
c     #        R = 1 + r*cos(2*pi*y)
c     # 
c     # The solution proceeds in two steps. If we are
c     # solving an incompressible problem, we only need
c     # to evolve (x(t), y(t)).  If we have a compressible
c     # problem, we first evolve (x(t),y(t)) back to a starting
c     # location (x0,y0), and then evolve x, y and q from 
c     # the starting location back to final position.
c     # --------------------------------------------------------

      double precision function qexact(blockno, x,y,tfinal)
      implicit none

      external swirl_rhs_divfree
      external swirl_rhs_nondivfree
      external solout

      double precision x,y,tfinal
      integer blockno

      double precision xc0, yc0
      double precision xc1, yc1

c      integer blockno_dummy
      double precision t0
      double precision xp,yp,zp

      double precision sigma(3), rtol, atol
      integer itol, iout

      integer Nmax, lwork,nrdens, liwork
      parameter(Nmax=3, nrdens=0)
      parameter(lwork=8*Nmax+5*nrdens+21)
      parameter(liwork=nrdens+21)

      double precision work(lwork), rpar
      double precision q0_physical, q0
      integer iwork(liwork), ipar(2), idid

      double precision tol

      logical evolve_q

      integer example
      common /example_comm/ example

      integer use_stream
      common /velocity_comm/ use_stream

      integer initchoice
      common /initchoice_comm/ initchoice


      integer*8 cont, get_context

      integer i

      cont = get_context()

c     # ------------------------------------------
c     # Numerical parameters
c     # ------------------------------------------
      itol = 0
      rtol = 1.d-12
      atol = 1.d-12
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

      ipar(1) = example
      ipar(2) = use_stream

      call fclaw_dopri5(2,swirl_rhs_divfree,t0,sigma,tfinal,
     &            rtol,atol,itol,
     &            solout,iout, work,lwork,iwork,liwork,
     &            rpar,ipar,idid)

      if (idid .ne. 1) then
          write(6,*) 'DOPRI5 : idid .ne. 1'
          stop
      endif

c     # Initial position traced back from (xc1,yc1)
      xc0 = sigma(1)
      yc0 = sigma(2)

c     # Map to [-1,1]x[-1,1]      
      call fclaw_map_2d_c2m_cart(blockno,xc0,yc0,xc1,yc1,zp)

c     # Map to [0,1]x[0,1]      
c      xp = (1 + xc0)/2.d0
c      yp = (1 + yc0)/2.d0
      xp = xc0
      yp = yc0

      if (initchoice .eq. 0) then
          write(6,*) 'qexact.f : initchoice .eq. 0 not defined'
          stop
      elseif (initchoice .eq. 1) then
          q0 = q0_physical(xp,yp,zp)
      elseif (initchoice .eq. 2) then
          q0 = 1.d0
      endif

      
      evolve_q = .true.
      if (evolve_q) then
c         # We now need to evolve q along with (x,y), starting from
c         # from (xc0,yc0)
          sigma(1) = xc0
          sigma(2) = yc0
          sigma(3) = q0

          do i = 1,20
              work(i) = 0
              iwork(i) = 0
          enddo

          t0 = 0
          call fclaw_dopri5(3,swirl_rhs_nondivfree,t0,sigma,tfinal,
     &                rtol,atol,itol,
     &                solout,iout, work,lwork,iwork,liwork,
     &                rpar,ipar,idid)

          if (idid .ne. 1) then
              write(6,*) 'DOPRI5 : idid .ne. 1'
              stop
          endif

          tol = 1e-8
          if (abs(sigma(1)-x) .gt. tol) then
              write(6,*) 'qexact.f : Did not evolve x correctly'
              write(6,100) xc0,yc0
              write(6,100) x,y
              write(6,100) sigma(1), sigma(2)
              write(6,105) abs(x-sigma(1)), abs(y-sigma(2))
              stop
          endif
          if (abs(sigma(2)-y) .gt. tol) then
              write(6,*) 'qexact.f : Did not evolve y correctly'
              write(6,100) xc0,yc0
              write(6,100) x,y
              write(6,100) sigma(1), sigma(2)
              write(6,105) abs(x-sigma(1)), abs(y-sigma(2))
              stop
          endif
100       format(2F24.16)
105       format(2E24.4)          

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


c     # ---------------------------------------------------------------      
      double precision function q0_physical(xp,yp,zp)
      implicit none

      double precision xp, yp, zp

      double precision x0, y0, z0, q0, r0, r
      double precision xpp, ypp
      double precision Hsmooth

c     # Sphere centered at (0.5,0.5,0) on swirl
      x0 = 0.5
      y0 = 0.5
      z0 = 0
      r0 = 0.2   !! radius of sphere

c     # modulo(A,P) = A - floor(A/P)*P
c      write(6,*) 'mod(6.1,1.d0)  ', modulo(6.1d0,1.d0)
c      write(6,*) 'mod(0.5,1.d0)  ', modulo(0.5d0,1.d0)
c      write(6,*) 'mod(-0.5,1.d0) ', modulo(-0.5d0,1.d0)
c      write(6,*) 'mod(-1.3,1.d0) ', modulo(-1.3d0,1.d0)
c      write(6,*) ' '
c      write(6,*) 'mod(6.1,1.d0)  ', 6.1d0 - floor(6.1d0)
c      write(6,*) 'mod(0.5,1.d0)  ', 0.5d0 - floor(0.5d0)
c      write(6,*) 'mod(-0.5,1.d0) ', -0.5d0 - floor(-0.5d0)
c      write(6,*) 'mod(-1.3,1.d0) ', -1.3d0 - floor(-1.3d0)
c      stop

      xpp = modulo(xp,1.d0)
      ypp = modulo(yp,1.d0)


      r = sqrt((xpp - x0)**2 + (ypp-y0)**2)
      q0 = Hsmooth(r + r0) - Hsmooth(r - r0)

      q0_physical = q0

      end

c     # ---------------------------------------------------------------      
      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end

