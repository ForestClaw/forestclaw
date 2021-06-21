c     
c     # Basic clawpack routine - with 4 additional arguments
c     # mwaves, mcapa, method, mthlim + ierror
c     ==========================================================
      subroutine clawpack46_step2(maxm,maxmx,maxmy,meqn,maux,mbc,
     &     mx,my,qold,aux,dx,dy,dt,cflgrid,fm,fp,gm,gp,
     &     faddm,faddp,gaddm,gaddp,q1d,dtdx1d,dtdy1d,
     &     aux1,aux2,aux3,work,mwork,rpn2,rpt2,flux2,
     &     mwaves,mcapa,method,mthlim, block_corner_count, ierror)
c     ==========================================================
      implicit none

      external rpn2, rpt2, flux2

      integer :: maxm, maxmx, maxmy, meqn, maux, mbc, mx, my, mwork
      double precision :: dx, dy, dt, cflgrid
      integer :: mwaves, mcapa, method(7), mthlim(mwaves)
      integer :: block_corner_count(0:3)
      integer :: ierror


      double precision :: qold(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision ::   fm(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision ::   fp(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision ::   gm(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision ::   gp(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision ::  q1d(1-mbc:maxm+mbc, meqn)
      double precision :: faddm(1-mbc:maxm+mbc, meqn)
      double precision :: faddp(1-mbc:maxm+mbc, meqn)
      double precision :: gaddm(1-mbc:maxm+mbc, meqn, 2)
      double precision :: gaddp(1-mbc:maxm+mbc, meqn, 2)
      double precision :: aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
      double precision :: aux1(1-mbc:maxm+mbc, maux)
      double precision :: aux2(1-mbc:maxm+mbc, maux)
      double precision :: aux3(1-mbc:maxm+mbc, maux)
      double precision :: dtdx1d(1-mbc:maxm+mbc)
      double precision :: dtdy1d(1-mbc:maxm+mbc)
      double precision :: work(mwork)

      double precision :: dtcom, dxcom, dycom, tcom
      integer :: icom, jcom
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      integer :: i0wave, i0s, i0amdq, i0apdq, i0cqxx, i0bmadq
      integer :: i0bpadq, iused
      double precision :: dtdx, dtdy, cfl1d
      integer :: m,i,j, ma, ixy
      integer :: sweep_dir

      integer :: blockno, fc2d_clawpack46_get_block

      ierror = 0
c     
c     # store mesh parameters that may be needed in Riemann solver but not
c     # passed in...
      dxcom = dx
      dycom = dy
      dtcom = dt

c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c     
      i0wave = 1
      i0s = i0wave + (maxm+2*mbc)*meqn*mwaves
      i0amdq = i0s + (maxm+2*mbc)*mwaves
      i0apdq = i0amdq + (maxm+2*mbc)*meqn
      i0cqxx = i0apdq + (maxm+2*mbc)*meqn
      i0bmadq = i0cqxx + (maxm+2*mbc)*meqn
      i0bpadq = i0bmadq + (maxm+2*mbc)*meqn
      iused = i0bpadq + (maxm+2*mbc)*meqn - 1
c     
      if (iused.gt.mwork) then
         ierror = 1
         return
      endif
c     
c     
      cflgrid = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
c     
      do  m=1,meqn
         do i=1-mbc,mx+mbc
            do j=1-mbc,my+mbc
               fm(i,j,m) = 0.d0
               fp(i,j,m) = 0.d0
               gm(i,j,m) = 0.d0
               gp(i,j,m) = 0.d0
            end do
         end do
      end do
c     
      if (mcapa.eq.0) then
c     # no capa array:
         do  i=1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
         end do
      endif

c     # Cubed sphere : Set corners for an x-sweep
c     # This does nothing for non-cubed-sphere grids. 
      sweep_dir = 0
      call clawpack46_fix_corners(mx,my,mbc,meqn,qold,sweep_dir,
     &     block_corner_count)


c     # perform x-sweeps
c     ==================
c     
      blockno = fc2d_clawpack46_get_block()
c      write(6,*) 'Doing x=sweep : blockno = ', blockno
      do  j = 2-mbc,my+mbc-1
c     
c     # copy data along a slice into 1d arrays:
         do m=1,meqn
            do i = 1-mbc, mx+mbc
               q1d(i,m) = qold(i,j,m)
            enddo
         enddo

         if (mcapa .gt. 0)  then
            do i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(i,j,mcapa)
            end do
         endif
c     
         if (maux .gt. 0)  then
            do  ma=1,maux
               do  i = 1-mbc, mx+mbc
                  aux1(i,ma) = aux(i,j-1,ma)
                  aux2(i,ma) = aux(i,j  ,ma)
                  aux3(i,ma) = aux(i,j+1,ma)
               end do
            end do
         endif
c     
c     
c     # Store the value of j along this slice in the common block
c     # comxyt in case it is needed in the Riemann solver (for
c     # variable coefficient problems)
         jcom = j
c     
c     # compute modifications fadd and gadd to fluxes along this slice:
         ixy = 1
         call flux2(1,maxm,meqn,maux,mbc,mx,
     &        q1d,dtdx1d,aux1,aux2,aux3,
     &        faddm,faddp,gaddm,gaddp,cfl1d,
     &        work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &        work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2,
     &        mwaves,mcapa,method,mthlim)

         cflgrid = dmax1(cflgrid,cfl1d)
c     
c     # update fluxes for use in AMR:
c     # NOTE : We update ghost cell values for subcycling
         do  m=1,meqn
            do  i=2-mbc,mx+mbc
               fm(i,j,m) = fm(i,j,m) + faddm(i,m)
               fp(i,j,m) = fp(i,j,m) + faddp(i,m)
               gm(i,j,m) = gm(i,j,m) + gaddm(i,m,1)
               gp(i,j,m) = gp(i,j,m) + gaddp(i,m,1)
               gm(i,j+1,m) = gm(i,j+1,m) + gaddm(i,m,2)
               gp(i,j+1,m) = gp(i,j+1,m) + gaddp(i,m,2)
            end do
         end do
      end do
c     
c     
c     
c     # perform y sweeps
c     ==================
c     

c     # Cubed sphere : Set corners for an y-sweep
c     # This does nothing for non-cubed-sphere grids. 
      sweep_dir = 1
      call clawpack46_fix_corners(mx,my,mbc,meqn,qold,sweep_dir,
     &     block_corner_count)

c     
      do  i = 2-mbc, mx+mbc-1
c     
c     # copy data along a slice into 1d arrays:
         do m=1,meqn
            do j = 1-mbc, my+mbc
               q1d(j,m) = qold(i,j,m)
            enddo
         enddo
c     
         if (mcapa .gt. 0)  then
            do  j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(i,j,mcapa)
            end do
         endif
c     
         if (maux .gt. 0)  then
            do  ma=1,maux
               do  j = 1-mbc, my+mbc
                  aux1(j,ma) = aux(i-1,j,ma)
                  aux2(j,ma) = aux(i,  j,ma)
                  aux3(j,ma) = aux(i+1,j,ma)
               end do
            end do
         endif
c     
c     
c     # Store the value of i along this slice in the common block
c     # comxyt in case it is needed in the Riemann solver (for
c     # variable coefficient problems)
         icom = i
c     
c     # compute modifications fadd and gadd to fluxes along this slice:
         ixy = 2
         call flux2(2,maxm,meqn,maux,mbc,my,
     &        q1d,dtdy1d,aux1,aux2,aux3,
     &        faddm,faddp,gaddm,gaddp,cfl1d,
     &        work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &        work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2,
     &        mwaves,mcapa,method,mthlim)

         cflgrid = dmax1(cflgrid,cfl1d)
c     
c     #
c     # update fluxes for use in AMR:
c     
         do  m=1,meqn
            do  j=2-mbc,my+mbc
               gm(i,j,m) = gm(i,j,m) + faddm(j,m)
               gp(i,j,m) = gp(i,j,m) + faddp(j,m)
               fm(i,j,m) = fm(i,j,m) + gaddm(j,m,1)
               fp(i,j,m) = fp(i,j,m) + gaddp(j,m,1)
               fm(i+1,j,m) = fm(i+1,j,m) + gaddm(j,m,2)
               fp(i+1,j,m) = fp(i+1,j,m) + gaddp(j,m,2)
            end do
         end do
      end do
c     
c     
      return
      end

c     #  See 'cubed_sphere_corners.ipynb'
      subroutine clawpack46_fix_corners(mx,my,mbc,meqn,q,sweep_dir,
     &     block_corner_count)
      implicit none

      integer :: mx,my,mbc,meqn,sweep_dir
      integer :: block_corner_count(0:3)
      double precision :: q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer :: k,m,idata,jdata
      double precision :: ihat(0:3),jhat(0:3)
      integer :: i1, j1, ibc, jbc
      logical :: use_b

c     # Lower left corner
      ihat(0) = 0.5
      jhat(0) = 0.5

c     # Lower right corner
      ihat(1) = mx+0.5
      jhat(1) = 0.5

c     # Upper left corner
      ihat(2) = 0.5
      jhat(2) = my+0.5

c     # Upper right corner
      ihat(3) = mx+0.5
      jhat(3) = my+0.5

      do k = 0,3
         if (block_corner_count(k) .ne. 3) then
            cycle
         endif
         use_b = sweep_dir .eq. 0 .and. (k .eq. 0 .or. k .eq. 3)
     &        .or.  sweep_dir .eq. 1 .and. (k .eq. 1 .or. k .eq. 2)
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Average fine grid corners onto coarse grid ghost corners
               if (k .eq. 0) then
                  i1 = 1-ibc
                  j1 = 1-jbc
               elseif (k .eq. 1) then
                  i1 = mx+ibc
                  j1 = 1-jbc
               elseif (k .eq. 2) then
                  i1 = 1-ibc
                  j1 = my+jbc
               elseif (k .eq. 3) then
                  i1 = mx+ibc
                  j1 = my+jbc
               endif

               if (use_b) then
c                 # Transform involves B                
                  idata =  j1 + int(ihat(k) - jhat(k))
                  jdata = -i1 + int(ihat(k) + jhat(k))
               else
c                 # Transform involves B.transpose()             
                  idata = -j1 + int(ihat(k) + jhat(k))
                  jdata =  i1 - int(ihat(k) - jhat(k))
               endif 
               do m = 1,meqn
                  q(i1,j1,m) = q(idata,jdata,m)
               end do           !! jbc
            end do              !! ibc
         end do                 !! corner 'k' loop
      end do                    !! meqn loop

      end
