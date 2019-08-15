      subroutine square_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer initchoice
      common /initchoice_comm/ initchoice

      integer example
      common /example_comm/ example      

      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      integer i,j, mq
      double precision qmin, qmax, xc, yc
      double precision xp,yp,zp

      logical refine
      


      tag_patch = 0

      cont = get_context()

c     # Refine based only on first variable in system.
      refine = .false.
      mq = 1
      qmin = q(1,1,mq)
      qmax = q(1,1,mq)
      do j = 1,my
         do i = 1,mx
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
            if (initchoice .le. 1) then
               if (example .eq. 0) then
                   refine = q(i,j,mq) .gt. tag_threshold .and. 
     &                q(i,j,mq) .lt. 1-tag_threshold
               else
                    refine = q(i,j,mq) .gt. tag_threshold
               endif
            elseif (initchoice .eq. 2) then
                refine = q(i,j,mq) .gt.  tag_threshold              
            else                
                write(6,'(A,A)') 'Refining not yet defined for ',
     &                  'example > 0'
                stop
            endif
            if (refine) then
                tag_patch = 1
                return
            endif
         enddo
      enddo

      end
