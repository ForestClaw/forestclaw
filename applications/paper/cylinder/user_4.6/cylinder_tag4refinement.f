      subroutine cylinder_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example      

      integer refine_pattern
      common /refine_comm/ refine_pattern

      integer initchoice
      common /initchoice_comm/ initchoice

      integer*8 cont, get_context

      integer i,j, mq
      double precision qmin, qmax, xc, yc
      double precision r, ravg, xp,yp,zp, theta, phi
      double precision xc1, yc1, zc1
      logical refine

      cont = get_context()

      if (blockno .ne. 0) then
          write(6,'(A,A)') 'cylinder_tag3refinement : Cylinder ',
     &           'tagging assumes single block'
          stop
      endif

      refine = .false.
      mq = 1
      do j = 1,my
         do i = 1,mx
            if (refine_pattern .eq. 0) then
                refine = (q(i,j,mq) .gt.  tag_threshold) .and.
     &                   (q(i,j,mq) .lt. 1-tag_threshold)
            else
                xc = xlower + (i-0.5)*dx
                yc = ylower + (j-0.5)*dy
                if (refine_pattern .eq. 1) then
                    refine = xc .gt. 0.5
                elseif (refine_pattern .eq. 2) then
                    refine = yc .gt. 0.5
                endif
            endif
            if (refine) then
                tag_patch = 1
                return
            endif
         enddo
      enddo


      end
