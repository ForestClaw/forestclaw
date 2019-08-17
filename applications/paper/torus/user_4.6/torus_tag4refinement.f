      subroutine torus_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer example
      common /example_comm/ example      

      integer refine_pattern
      common /refine_comm/ refine_pattern

      integer initchoice
      common /initchoice_comm/ initchoice

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      integer*8 cont, get_context

      integer i,j, mq
      double precision qmin, qmax, xc, yc
      double precision r, ravg, xp,yp,zp
      logical refine

      cont = get_context()

      refine = .false.
      mq = 1
      do j = 1,my
         do i = 1,mx
            if (refine_pattern .eq. 0) then
                refine = q(i,j,mq) .gt.  tag_threshold              
            else
                xc = xlower + (i-0.5)*dx
                yc = ylower + (j-0.5)*dy
                call fclaw2d_map_c2m(cont,blockno,xc,yc,
     &                                      xp, yp,zp);
                if (refine_pattern .eq. 1) then
                    refine = xp .lt. 0
                elseif (refine_pattern .eq. 2) then
                    r = sqrt(xp**2 + yp**2)
                    refine = r .gt. 1.d0
                endif
            endif
            if (refine) then
                tag_patch = 1
                return
            endif
         enddo
      enddo


      end

