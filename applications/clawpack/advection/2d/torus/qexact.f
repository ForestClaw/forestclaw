c     # ---------------------------------------------------------------      
      double precision function qexact(blockno,xc,yc,t)
      implicit none

      external psi_rhs_torus
      external solout

      integer blockno
      double precision xc,yc,t

      integer*8 cont, get_context, blockno_dummy
      double precision xc1, yc1, zc1, t0, tfinal
      double precision xp,yp,zp

      double precision sigma(2), rtol, atol
      integer itol, iout

      integer Nmax, lwork,nrdens, liwork
      parameter(Nmax=2, nrdens=0)
      parameter(lwork=8*Nmax+5*nrdens+21)
      parameter(liwork=nrdens+21)

      double precision work(lwork), rpar
      double precision q0_physical, q0
      integer iwork(liwork), ipar, idid

      integer i

      sigma(1) = xc
      sigma(2) = yc

      cont = get_context()

      itol = 0
      rtol = 1.d-4
      atol = 1.d-4
      iout = 0

      do i = 1,20
          work(i) = 0
          iwork(i) = 0
      enddo


c     # This is not the torus mapping, but rather maps the brick to
c     # a unit square      
      call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

      sigma(1) = xc1
      sigma(2) = yc1
      t0 = 0
      tfinal = t
      call dopri5(2,psi_rhs_torus,t0,sigma,tfinal,rtol,atol,itol,
     &            solout,iout, work,lwork,iwork,liwork,
     &            rpar,ipar,idid)

      if (idid .ne. 1) then
          write(6,*) 'DOPRI5 : idid .ne. 1'
          stop
      endif

      xc1 = sigma(1)
      yc1 = sigma(2)

      !! Do not map (xc1,yc1) to brick, since mapping was done above
      blockno_dummy = -1  

      call fclaw2d_map_c2m(cont,blockno_dummy,xc1,yc1,xp,yp,zp)

      q0 = q0_physical(xp,yp,zp)

      qexact = q0

      end

c     # ---------------------------------------------------------------      
      subroutine solout(nr,xold,x,y,n,con,icomp,nd,rpar,ipar,irtrn)
      dimension y(n),con(5*nd),icomp(nd)

c     # Dummy routine

      end


c     # ---------------------------------------------------------------      
      double precision function q0(blockno,xc1,yc1)
      implicit none

      integer blockno
      double precision xc1, yc1

      integer*8 cont, get_context
      double precision xp, yp, zp
      double precision q0_physical

      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc1,yc1,xp,yp,zp)

      q0 = q0_physical(xp,yp,zp)

      end


c     # ---------------------------------------------------------------      
      double precision function q0_physical(xp,yp,zp)
      implicit none

      double precision xp, yp, zp

      double precision r,r0, x0, y0, z0, q0
      double precision Hsmooth

      integer mapping
      common /mapping_comm/ mapping

      double precision alpha, revs_per_s
      common /torus_comm/ alpha, revs_per_s      

c     # Sphere centered at (1,0,r0) on torus
      r0 = alpha
      if (mapping .le. 1) then
          x0 = 1.0
          y0 = 0.0
          z0 = r0
      else
          x0 = 0.5
          y0 = 0.5
          z0 = 0
      endif

      r = sqrt((xp - x0)**2 + (yp-y0)**2 + (zp-z0)**2)
      q0 = Hsmooth(r + r0) - Hsmooth(r - r0)

      q0_physical = q0

      end

c     # ---------------------------------------------------------------      
      double precision function Hsmooth(r)
      implicit none

      double precision r

      Hsmooth = (tanh(r/0.02d0) + 1)/2.d0

      end
