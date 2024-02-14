c--------------------------------------------------------------------
c> modified interpolation stencil that grabs from corners
c--------------------------------------------------------------------
      subroutine test_2d_clawpatch46_fort_interpolate_face
     &      (mx,my,mbc,meqn,qcoarse,qfine,
     &      idir,iface_coarse,num_neighbors,refratio,igrid,
     &      transform_ptr)
      implicit none
      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer num_neighbors
      integer*8 transform_ptr
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq,r2, m
      integer ibc, i1
      integer jbc, j1
      integer ic, jc, mth
      double precision qc, qo

c     # This should be refratio*refratio.
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      logical fclaw2d_clawpatch_is_valid_interp
      logical skip_this_grid

      integer a(2,2), f(2)
      integer ii,jj,dc(2),df(2,0:rr2-1),iff,jff
      integer shiftx(0:rr2-1),shifty(0:rr2-1)

      mth = 5
      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif


      call fclaw2d_clawpatch_build_transform(transform_ptr,a,f)

c     # This needs to be written for refratios .ne. 2.
      m = 0
      do jj = 0,1
         do ii = 0,1
c           # Direction on coarse grid
            dc(1) = ii
            dc(2) = jj

c           # Direction on fine grid (converted using metric). Divide
c           # by refratio to scale length to unit vector
            df(1,m) = (a(1,1)*dc(1) + a(1,2)*dc(2))/refratio
            df(2,m) = (a(2,1)*dc(1) + a(2,2)*dc(2))/refratio

c           # Map (0,1) to (-1/4,1/4) (locations of fine grid points)
            shiftx(m) = -1 + 2*ii
            shifty(m) = -1 + 2*jj
            m = m + 1
         enddo
      enddo
c     # Create map :

      do mq = 1,meqn
         if (idir .eq. 0) then
c           # this ensures that we get 'hanging' corners

            do ibc = 1,mbc/2
            if (iface_coarse .eq. 0) then
               ic = ibc
            elseif (iface_coarse .eq. 1) then
               ic = mx - ibc + 1
            else
               write(6,*) 'interpolate : Problem with iface_coarse'
               write(6,*) 'iface_coarse = ', iface_coarse
               stop               
            endif
            do jc = 1,mx
               i1 = ic
               j1 = jc
               call fclaw2d_clawpatch_transform_face_half(i1,j1,i2,j2,
     &               transform_ptr)
               skip_this_grid = .false.
               do m = 0,r2-1
                  if (.not. 
     &        fclaw2d_clawpatch_is_valid_interp(i2(m),j2(m),mx,my,mbc))
     &                  then
                     skip_this_grid = .true.
                     exit
                  endif
               enddo
               if (.not. skip_this_grid) then
                  qc = qcoarse(ic,jc,mq)
                  do m = 0,rr2-1
                     qo = qcoarse(ic+shiftx(m),jc+shifty(m),mq)
                     iff = i2(0) + df(1,m)
                     jff = j2(0) + df(2,m)
                     qfine(iff,jff,mq) = 0.75d0*qc + 0.25d0*qo
                  enddo
               endif
            enddo
            enddo
         else
            do jbc = 1,mbc/2
            if (iface_coarse .eq. 2) then
               jc = jbc
            elseif (iface_coarse .eq. 3) then
c              // iface_coarse = 3
               jc = my - jbc + 1
            else
               write(6,*) 'interpolate : Problem with iface_coarse'
               write(6,*) 'iface_coarse = ', iface_coarse
               stop
            endif
            do ic = 1,mx
               i1 = ic
               j1 = jc
               call fclaw2d_clawpatch_transform_face_half(i1,j1,i2,j2,
     &               transform_ptr)
c              # ---------------------------------------------
c              # Two 'half-size' neighbors will be passed into
c              # this routine.  Only half of the coarse grid ghost
c              # indices will be valid for the particular grid
c              # passed in.  We skip those ghost cells that will
c              # have to be filled in by the other half-size
c              # grid.
c              # ---------------------------------------------
               skip_this_grid = .false.
               do m = 0,r2-1
                  if (.not. 
     &       fclaw2d_clawpatch_is_valid_interp(i2(m),j2(m),mx,my,mbc))
     &                  then
                     skip_this_grid = .true.
                     exit
                  endif
               enddo
               if (.not. skip_this_grid) then
                  qc = qcoarse(ic,jc,mq)
                  do m = 0,rr2-1
                     qo = qcoarse(ic+shiftx(m),jc+shifty(m),mq)
                     iff = i2(0) + df(1,m)
                     jff = j2(0) + df(2,m)
                     qfine(iff,jff,mq) = 0.75d0*qc + 0.25d0*qo
                  enddo

               endif                    !! Don't skip this grid
            enddo                       !! i loop
            enddo                       !! end of jbc loop
         endif                          !! end idir branch
      enddo                             !! endo mq loop

      end

c--------------------------------------------------------------------
c> modified interpolation stencil that grabs from corners
c--------------------------------------------------------------------
      subroutine test_2d_clawpatch46_fort_interpolate_corner
     &      (mx,my,mbc,meqn,refratio,
     &      qcoarse,qfine,icorner_coarse,transform_ptr)
      implicit none

      integer mx,my,mbc,meqn,icorner_coarse,refratio
      integer*8 transform_ptr
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ic, jc, mq, ibc,jbc, mth
      double precision qc, qo

c     # This should be refratio*refratio.
      integer i1,j1,m, r2
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)

      integer a(2,2), f(2)
      integer ii,jj,iff,jff,dc(2),df(2,0:rr2-1)
      integer shiftx(0:rr2-1), shifty(0:rr2-1)

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_corner_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif

      call fclaw2d_clawpatch_build_transform(transform_ptr,a,f)

      m = 0
      do jj = 0,1
         do ii = 0,1
c           # Direction on coarse grid
            dc(1) = ii
            dc(2) = jj

c           # Direction on fine grid (converted using metric). Divide
c           # by 2 (refratio) to scale length to unit vector
            df(1,m) = (a(1,1)*dc(1) + a(1,2)*dc(2))/2
            df(2,m) = (a(2,1)*dc(1) + a(2,2)*dc(2))/2

c           # Map (0,1) to (-1/4,1/4) (locations of fine grid points)
            shiftx(m) = -1 + 2*ii
            shifty(m) = -1 + 2*jj
            m = m + 1
         enddo
      enddo


      mth = 5

      do ibc = 1,mbc/2
      do jbc = 1,mbc/2
      if (icorner_coarse .eq. 0) then
         ic = ibc
         jc = jbc
      elseif (icorner_coarse .eq. 1) then
         ic = mx - ibc + 1
         jc = jbc
      elseif (icorner_coarse .eq. 2) then
         ic = ibc
         jc = my - jbc + 1
      elseif (icorner_coarse .eq. 3) then
         ic = mx - ibc + 1
         jc = my - jbc + 1
      else
         write(6,*) "interpolate : Problem with icorner_coarse"
         write(6,*) "icorner_coarse = ", icorner_coarse
         stop
      endif

c     # Interpolate coarse grid corners to fine grid corner ghost cells
      i1 = ic
      j1 = jc
      call fclaw2d_clawpatch_transform_corner_half(i1,j1,i2,j2,
     &      transform_ptr)

      do mq = 1,meqn
         qc = qcoarse(ic,jc,mq)
         do m = 0,rr2-1
            qo = qcoarse(ic+shiftx(m),jc+shifty(m),mq)
            iff = i2(0) + df(1,m)
            jff = j2(0) + df(2,m)
            qfine(iff,jff,mq) = 0.75d0*qc + 0.25d0*qo
         enddo

      enddo
      enddo 
      enddo

      end


