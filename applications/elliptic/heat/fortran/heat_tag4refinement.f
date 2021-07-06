      subroutine tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qvec(3), dmax, heat_eval_refinement

      tag_patch = 0
c     # Refine based only on first variable in system.
      mq = 1
      do j = 1,my+1
         do i = 1,mx+1
            qvec(1) = q(i,j,1)
            qvec(2) = q(i-1,j,1)
            qvec(3) = q(i,j-1,1)

            dmax = heat_eval_refinement(qvec,dx,dy)

            if (abs(dmax) .gt. tag_threshold) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end      

      double precision function heat_eval_refinement(q,dx,dy)
          implicit none

          double precision q(3), dx,dy

          double precision qij, qimj, qijm
          double precision qxm, qym, dmax

          qij = q(1)
          qimj = q(2)
          qijm = q(3)
          qxm = (qij - qimj)/dx
          qym = (qij - qijm)/dy

          dmax = max(abs(qxm),abs(qym))

          heat_eval_refinement = dmax

          return

      end function heat_eval_refinement

