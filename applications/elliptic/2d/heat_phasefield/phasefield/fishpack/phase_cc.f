      program phase
      implicit none


c     # Grid parameters
      integer maxmx, maxmy
      parameter(maxmx = 2**12, maxmy = 2**12)


c     # Model parameters
      double precision xi, Tm, Lv, C, kappa, gma, mu, k,
     &      sigma, tau, beta, S, alpha, m, gma2, beta_p,
     &      tau_th, tanh

c     # Temporary variables
      double precision uij_n, src, d, dpdt, gphi,
     &  pij_n, pij_np1

c     # Input parameters
      double precision ax_in, ay_in, dx_in, dy_in, t_in, tstart
      double precision domain_length
      integer mx_in, my_in, nstp, izero, nchar
      character*100 fname
      logical restart

      double precision phi(maxmx, maxmy)
      double precision u(maxmx, maxmy)
      double precision theta, r, rp, s1, s2, amp

c     # More grid parameters
      integer mx, my, i, j, n, nout, nstep
      double precision ax, ay, bx, by, dx, dy, dt, tfinal, t
      double precision x(maxmx), y(maxmy)

c     # Misc
      integer iframe, iframe_start
      double precision pi, dsgn

      integer mbdcnd, nbdcnd, ierror, idimf
      double precision lambda_u, lambda_phi, pertrb
      double precision bda(maxmy), bdb(maxmy)
      double precision bdc(maxmx), bdd(maxmx)
      double precision f(maxmx, maxmy)

      double precision r0, x0, y0, kanio, pdel
      common /parm_comm/ r0, x0, y0, kanio, pdel

      integer wsize, wsize_t
c         13M + 4N + M*INT(LOG2(N))
c             log2(N) <= log2(maxmy) ~ 10
      parameter(wsize = 4*maxmy + (13 + 12)*maxmx)
      double precision work(wsize)

      double precision xlow, ylow, w


      wsize_t = 4*(maxmy) +
     &      (13 + int(log(real(my))/log(2.d0))*mx)
      if (wsize < wsize_t) then
         write(6,*) 'Increase wsize from ', wsize, ' to ',wsize_t
         stop
      endif
      pi = 4.d0*datan(1.d0)


c     # --------------------------------
c     # Set up grid
c     # --------------------------------

      open(10,file='phase.dat')
      read(10,*) mx, my
c      read(10,*) ax, bx, ay, by
      read(10,*) domain_length
      read(10,*) nstep, nout
      read(10,*) dt

      if (mx > maxmx .or. my > maxmy) then
         write(6,*) 'mx or my too big'
         stop
      endif

      ax = -domain_length/2
      bx = domain_length/2
      ay = -domain_length/2
      by = domain_length/2


      dx = (bx - ax)/mx
      dy = (by - ay)/my

      do i = 1,mx
         x(i) = ax + (i-0.5)*dx
      enddo

      do j = 1,my
         y(j) = ay + (j-0.5)*dy
      enddo

c     # ----------------------------
c     # physical parameters
c     # ----------------------------
      read(10,*) tm
      read(10,*) lv
      read(10,*) C
      read(10,*) kappa
      read(10,*) gma
      read(10,*) mu
      read(10,*) k
      read(10,*) sigma
      read(10,*) tau
      read(10,*) beta

c     # ----------------------------
c     # model parameters
c     # ----------------------------
      read(10,*) s
      read(10,*) alpha
      read(10,*) m
      read(10,*) xi
      read(10,*) kanio
      read(10,*) r0
      read(10,*) amp
      read(10,*) pdel
      read(10,*) gma2
      read(10,*) beta_p
      read(10,*) tau_th

c     # ----------------------------
c     # Get restart information
c     # ----------------------------
      read(10,*) restart
      if (restart) then
         read(10,*) fname
      endif
      close(10)

      write(6,*) 'Parameters used for this run : '
      write(6,*) '-------------------------------'
      write(6,*) 'Physical parameters '
      write(6,*) '-------------------------------'
      write(6,200) 'Tm    = ', tm
      write(6,200) 'Lv    = ', Lv
      write(6,200) 'C     = ', C
      write(6,200) 'kappa = ', kappa
      write(6,210) 'Gamma = ', gma
      write(6,200) 'mu    = ', mu
      write(6,200) 'k     = ', k
      write(6,210) 'sigma = ', sigma
      write(6,210) 'tau   = ', tau
      write(6,200) 'beta  = ', beta
      write(6,*) ' '
      write(6,*) '---------------------------'
      write(6,*) 'Non-dimensional parameters'
      write(6,*) '---------------------------'
      write(6,200) 'S       = ', S
      write(6,200) 'alpha   = ', alpha
      write(6,200) 'm       = ', m
      write(6,200) 'xi      = ', xi
      write(6,200) 'k       = ', kanio
      write(6,200) 'r0      = ', r0
      write(6,200) 'amp     = ', amp
      write(6,200) 'pdel    = ', pdel
      write(6,200) 'Gamma   = ', gma2
      write(6,200) 'beta_p  = ', beta_p
      write(6,200) 'tau_th  = ', tau_th
      write(6,*) ' '
      write(6,*) '---------------------------'
      write(6,*) 'Grid parameters'
      write(6,*) '---------------------------'
      write(6,220) 'mx   = ', mx
      write(6,220) 'my   = ', my
      write(6,200) 'ax   = ', ax
      write(6,200) 'bx   = ', bx
      write(6,200) 'ay   = ', ay
      write(6,200) 'by   = ', by
      write(6,200) 'dx   = ', dx
      write(6,200) 'dy   = ', dy

      if (restart) then
         write(6,*) 'Restart from file ', fname
      else
         write(6,*) 'Initializing run from scratch...'
      endif

 200  format(a,F10.5)
 210  format(a,E12.6)
 220  format(a,I5)

c     # --------------------------------
c     # Initialize phi
c     # --------------------------------
      write(6,*) 'Initializing phi and u...'

      if (restart) then
         open(20,file=fname)
         read(20,*) mx_in
         read(20,*) my_in
         read(20,*) ax_in
         read(20,*) ay_in
         read(20,*) dx_in
         read(20,*) dy_in
         read(20,*) t_in

         if (mx_in .ne. mx .or. my_in .ne. my) then
            write(6,*) 'Problem with restart: mx or my not right.'
            stop
         elseif (ax_in .ne. ax .or. ay_in .ne. ay) then
            write(6,*) 'Problem with restart: ax or ay not right'
            stop
         endif
         izero = iachar('0')
         tstart = t_in
         iframe_start = 0
         nstp = 1
         do i = 10,7,-1
            nchar = iachar(fname(i:i)) - izero
            iframe_start = iframe_start + nchar*nstp
            nstp = nstp*10
         enddo

         do j = 1,my
            do i = 1,mx
               read(20,*) phi(i,j), u(i,j)
            enddo
         enddo
         close(20)
      else

         x0 = (ax + bx)/2.d0
         y0 = (ay + by)/2.d0
         tstart = 0.0
         iframe = 0
         do j = 1,my
            do i = 1,mx
                xlow = ax + (i-1)*dx
                ylow = ay + (j-1)*dy
c                rp = sqrt((xlow)**2 + (ylow)**2)
 
                call cellave(xlow,ylow,dx,dy,w)
                u(i,j) = w-1.d0
                phi(i,j) = 1-w

c                r = sqrt((amp*(x(i) - x0))**2 + (y(j) - y0)**2)
c               theta = atan2(y(j) - y0,x(i) - x0)
c               rp = r0*(1.d0 + pdel*cos(kanio*theta))
c               d = tanh((r-rp)/xi)
c               u(i,j) = -0.5d0*(d + 1)  !! In [-1,0]
c               phi(i,j) = 0 0.5d0*(d + 1)  !! In [0,1]
            enddo
         enddo
         write(6,*) 'Done with initialization'

c     # Output frame only if we have initialized data from scratch.
         write(6,600) 0,tstart
         call out2(maxmx,maxmy,mx,my,ax,ay,dx,dy,phi,u,tstart,
     &        iframe_start)
      endif


c     % --------------------------------------
c     % Time step through results
c     % --------------------------------------

      do j = 1,my
         bda(j) = 0.d0
         bdb(j) = 0.d0
      enddo

      do i = 1,mx
         bdc(i) = 0.d0
         bdd(i) = 0.d0
      enddo

      lambda_phi = -1.d0/(m*dt)
      lambda_u = -1.d0/dt
      idimf = maxmx
      mbdcnd = 3
      nbdcnd = 3

      iframe = iframe_start
      do n = 1,nout
         t = tstart + n*dt

c         write(6,100) n,t
  100    format('Step',I5,' at time t = ',1PE12.5)

c        # Compute right hand side for phi :
         do j = 1,my
            do i = 1,mx
               pij_n = phi(i,j)
               uij_n = u(i,j)
               gphi = pij_n**2*(1.d0-pij_n)**2
               src = 30.d0*xi*S*alpha*gphi*uij_n +
     &              pij_n*(1.d0 - pij_n)*(pij_n - 0.5d0)
               f(i,j) = lambda_phi*(dt*(m/(xi**2))*src + pij_n)
            enddo
         enddo

c        # Solve equation for phi; solution will be in f
         call hstcrt(ax,bx,mx,mbdcnd,bda,bdb,ay,by,my,nbdcnd,
     &         bdc,bdd,lambda_phi,f,idimf,pertrb,ierror,work)

         if (ierror .ne. 0) then
            write(6,*) 'hstcrt : ierror /= 0; ', ierror
            stop
         endif

         if (work(1) > wsize) then
            write(6,*) 'hstcrt : work(1) > wsize; ', work(1), wsize
            stop
         endif

         if (abs(pertrb) .gt. 0) then
            write(6,*) 'pertrb > 0;  ', pertrb
            stop
         endif

c        # Set up right hand side for u
         do j = 1,my
            do i = 1,mx
               pij_n = phi(i,j) !! old phi(i,j)
               pij_np1 = f(i,j) !! new phi(i,j)
               phi(i,j) = f(i,j)  !! for next time step

               uij_n = u(i,j)
               dpdt = (pij_np1 - pij_n)/dt !! d(phi)/dt
               gphi = pij_np1**2*(1.0 - pij_np1)**2
               src = 30.d0*gphi*(-1.d0/S)*dpdt
               f(i,j) = lambda_u*(dt*src + uij_n)
            enddo
         enddo

c        # Solve equation for u; solution will be in f
         call hstcrt(ax,bx,mx,mbdcnd,bda,bdb,ay,by,my,nbdcnd,
     &        bdc,bdd,lambda_u,f,idimf,pertrb,ierror,work)

         if (ierror .ne. 0) then
            write(6,*) 'hstcrt : ierror /= 0', ierror
            stop
         endif

         if (abs(pertrb) .gt. 0) then
            write(6,*) 'pertrb > 0;  ', pertrb
            stop
         endif

         do j = 1,my
            do i = 1,mx
               u(i,j) = f(i,j)
            enddo
         enddo

c        # ---------------------------------
c        # Output results
c        # ---------------------------------

          write(6,601) n, t

         if (nstep*(n/nstep) == n) then
            iframe = iframe + 1
            write(6,600) iframe, t
            write(6,*) ' '
            call out2(maxmx,maxmy,mx,my,ax,ay,dx,dy,phi,u,t,iframe)
         endif
  600    format('Writing frame ',I5,' at time t = ',1PE12.5)
  601    format('Step ',I5,' at time t = ',1PE12.5)


      enddo

      end

      double precision function fdisc(x,y)
      implicit none

      double precision x,y, r, rp, theta

      double precision r0, x0, y0, kanio, pdel
      common /parm_comm/ r0, x0, y0, kanio, pdel

      r = sqrt((x-x0)**2 + (y-y0)**2)
      theta = atan2(y - y0,x - x0)
      rp = r0*(1.d0 + pdel*cos(kanio*theta))
      fdisc = r - rp

      end



c      double precision function tanh(x)
c      implicit none
c      double precision x
c      tanh = (exp(x) - exp(-x))/(exp(x) + exp(-x))
c      return
c      end
