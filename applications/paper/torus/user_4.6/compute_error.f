      subroutine torus46_compute_error(blockno, mx,my,mbc,meqn,
     &      dx,dy,xlower,ylower,t,q,error, soln)
      implicit none

      integer mx,my,mbc,meqn, blockno
      double precision dx, dy, xlower, ylower, t
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision soln(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer example
      common /example_comm/ example    

      integer*8 cont, get_context

      integer i,j,m
      double precision xc,yc, qexact
      double precision xc1, yc1, zc1, x,y

      double precision torus_divergence
      integer flow_flag
      
      cont = get_context()

      if (example .eq. 0) then
c          # non-divergent flow         
           flow_flag = 0
      else
c          # Divergent flow         
           flow_flag = 1
      endif

c     # Assume a single field variable only
      do j = 1,my
         do i = 1,mx
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            if (t .eq. 0) then
               soln(i,j,1) = q(i,j,1)
            else
               call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)
               soln(i,j,1) = qexact(xc1,yc1,t,flow_flag)
            endif
            error(i,j,1) = q(i,j,1) - soln(i,j,1)
         enddo
      enddo


      end
