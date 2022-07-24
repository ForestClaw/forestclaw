!!> @file
!!> Fortran subroutines for metric terms
!!
!! ---------------------------------------------------------
!!> @brief @copybrief ::fclaw3d_fort_compute_mesh_t
!!>
!!> Default implimentation
!!>
!!> @details @copydetails ::fclaw3d_fort_compute_mesh_t
!! ---------------------------------------------------------
subroutine fclaw3d_metric_fort_compute_mesh(mx,my,mz,mbc, &
        xlower,ylower, zlower, dx,dy,dz, blockno, & 
        xp,yp,zp,xd,yd,zd)
    implicit none

    integer mx,my, mz, mbc, blockno
    double precision dx,dy,dz, xlower,ylower, zlower

    double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

    double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)

    integer i,j,k
    double precision dxf,dyf, dzf, xc,yc,zc, xd1, yd1,zd1

    integer*8 map_context_ptr, fclaw_map_get_context

    map_context_ptr = fclaw_map_get_context()

    !! # We need both cell centered and node locations to
    !! # compute the basis at each face
    dxf = dx/2.d0
    dyf = dy/2.d0
    dzf = dz/2.d0
    do i = -2*mbc-1,2*(mx+mbc+1)+1
        do j = -2*mbc-1,2*(my+mbc+1)+1
            do k = -2*mbc-1,2*(mz+mbc+1) + 1
                !! Skip over any values at face centers
                if ((abs(mod(i,2)) .ne. abs(mod(j,2))) .or.  & 
                    (abs(mod(i,2)) .ne. abs(mod(k,2))) .or.  & 
                    (abs(mod(j,2)) .ne. abs(mod(k,2)))) then 
                    cycle
                endif
                !! Compute mesh value on 2x finer mesh; Choose some values for
                !! cell center on coarser mesh; some values for nodes.
                xc = xlower + (i-1)*dxf
                yc = ylower + (j-1)*dyf
                zc = zlower + (k-1)*dzf

                call fclaw3d_map_c2m(map_context_ptr, &
                                     blockno,xc,yc,zc,xd1,yd1,zd1)

                write(6,*) xc,yc,zc
                write(6,*) xd1, yd1, zd1
                write(6,*) ' '

                !! whole integer indices are cell centers. 
                if (abs(mod(i,2)) .eq. 1) then
                    !! # For odd values mesh values
                    !! # Physical location of cell vertices
                    xd((i-1)/2 + 1, (j-1)/2 + 1, (k-1)/2 + 1) = xd1
                    yd((i-1)/2 + 1, (j-1)/2 + 1, (k-1)/2 + 1) = yd1
                    zd((i-1)/2 + 1, (j-1)/2 + 1, (k-1)/2 + 1) = zd1
                else
                    !! # Physical locations of cell centers
                    xp(i/2,j/2,k/2) = xd1
                    yp(i/2,j/2,k/2) = yd1
                    zp(i/2,j/2,k/2) = zd1
                endif
            end do
        end do
    end do
    stop
end subroutine fclaw3d_metric_fort_compute_mesh

!! ----------------------------------------------------------
!!> @brief @copybrief ::fclaw3d_fort_compute_area_t
!!>
!!> Default implementation
!!>
!!> @details @copydetails ::fclaw3d_fort_compute_area_t
!! ---------------------------------------------------------
subroutine fclaw3d_metric_fort_compute_volume(mx,my,mz, mbc,dx,dy,dz, &
           xlower, ylower, zlower, blockno, xd,yd,zd, & 
           volume, faceareas, hexsize, hexstore, ghost_only)
    implicit none

    integer mx,my,mz,mbc,blockno, hexsize
    double precision hexstore(0:hexsize,0:hexsize,0:hexsize,3)
    double precision dx,dy, dz, xlower, ylower, zlower
    integer ghost_only
    double precision volume(-mbc:mx+mbc+1,-mbc:my+mbc+1, -mbc:mz+mbc+1)
    double precision faceareas(-mbc:mx+mbc+1,-mbc:my+mbc+1, -mbc:mz+mbc+2,3)

    double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)

    logical isaffine


    if (isaffine()) then
        !! # We don't need to compute areas all the way to the
        !! # finest level.
        call fclaw3d_metric_fort_compute_volume_affine(mx, my, mz, mbc, & 
                    xd, yd, zd, volume, faceareas, ghost_only)
    else
        call fclaw3d_metric_fort_compute_volume_general(mx, my, mz, mbc, & 
                dx, dy, dz, &
                xlower, ylower, zlower, blockno, volume,faceareas, &
                hexsize, hexstore, ghost_only)
    endif
end subroutine fclaw3d_metric_fort_compute_volume

!! ---------------------------------------------------------
!!> @brief Compute the volume for each cell for general mappings
!!>
!!> @param[in] mx, my, mz the number of cells in the x, y, z directions
!!> @param[in] mbc the number of ghost cells
!!> @param[in] dx, dy, dz the spacings in the x, y and z direcitons
!!> @param[in] xlower, ylower, zlower the lower left coordinate of the patch
!!> @param[in] blockno the block number
!!> @param[out] volume the volume of each cell
!!> @param[in] hexsize the length of the hex
!!> @param[in] hexstore stores a group of cell values
!!> @param[in] ghost_only
!! ---------------------------------------------------------
subroutine fclaw3d_metric_fort_compute_volume_general(mx,my,mz, mbc, & 
           dx, dy, dz, xlower, ylower, zlower, blockno, volume, &
           faceareas, hexsize, hexfine,ghost_only)
    implicit none

    integer :: mx,my,mz, mbc,blockno, hexsize
    double precision :: dx,dy, dz, xlower, ylower, zlower
    double precision :: hexfine(0:hexsize,0:hexsize,0:hexsize,3)
    integer :: ghost_only

    double precision volume(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    double precision faceareas(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3)

    integer i,j,k,ii,jj, kk, m, rfactor, icell, jcell, kcell
    double precision sum_volume
    double precision xp1,yp1,zp1
    double precision hex(0:1,0:1,0:1,3)
    double precision area(3), sum_face_area(3)

    double precision hex_compute_volume

    double precision dxf, dyf, dzf
    double precision xef, yef, zef, xe,ye, ze
    logical hex_is_volume_interior

    integer*8 cont, fclaw_map_get_context

    cont = fclaw_map_get_context()

    rfactor = hexsize
    dxf = dx/rfactor
    dyf = dy/rfactor
    dzf = dz/rfactor

    do j = -mbc,my+mbc+1
        do i = -mbc,mx+mbc+1
            do k = -mbc,mz+mbc+1

                !! Skip interior cells if we are only trying to rebuild
                !! metric terms in ghost cells
                if (hex_is_volume_interior(mx,my,mz,i,j,k) .and. & 
                    ghost_only .eq. 1) then
                    cycle
                endif
                xe = xlower + (i-1)*dx
                ye = ylower + (j-1)*dy
                ze = zlower + (k-1)*dz

                !! For each coarse grid cell, construct a local finer mesh at
                !! maximum level of refinement
                do ii = 0,rfactor
                    do jj = 0,rfactor
                        do kk = 0,rfactor
                            xef = xe + ii*dxf
                            yef = ye + jj*dyf
                            zef = ze + kk*dzf

                            call fclaw3d_map_c2m(cont, & 
                                blockno,xef,yef,zef, xp1,yp1,zp1)

                            hexfine(ii,jj,kk,1) = xp1
                            hexfine(ii,jj,kk,2) = yp1
                            hexfine(ii,jj,kk,3) = zp1
                        end do
                    end do
                end do

                !! Compute volume of each of the "sub-hexes"
                sum_volume = 0.d0
                do m = 1,3
                    sum_face_area(m) = 0.d0
                enddo
                do ii = 0,rfactor-1
                    do jj = 0,rfactor-1
                        do kk = 0, rfactor-1

                            !! Create sub-hex 
                            do icell = 0,1
                                do jcell = 0,1
                                    do kcell = 0,1
                                        do m = 1,3
                                            hex(icell,jcell,kcell,m) =  & 
                                              hexfine(ii+icell,jj+jcell,kk+kcell,m)
                                        end do
                                    end do
                                end do
                            end do

                            !! Compute volume of sub-hex at finest level
                            sum_volume = sum_volume + hex_compute_volume(hex)

                            !! Compute face areas by summing areas over sub-hexes
                            !! at faces which intersect coarser grid face
                            if (ii .eq. 0 .or. & 
                                jj .eq. 0 .or. & 
                                kk .eq. 0) then
                                call hex_compute_surf_area(hex,area)
                            endif

                            if (ii .eq. 0) then
                                sum_face_area(1) = sum_face_area(1) + area(1)
                            endif
                            if (jj .eq. 0) then
                                sum_face_area(2) = sum_face_area(2) + area(2)
                            endif
                            if (kk .eq. 0) then
                                sum_face_area(3) = sum_face_area(3) + area(3)
                            endif
                        end do
                    end do
                end do
                volume(i,j,k) = sum_volume

                do m = 1,3
                    faceareas(i,j,k,m) = sum_face_area(m)
                end do

            end do
        end do
    end do
end subroutine fclaw3d_metric_fort_compute_volume_general


!! ---------------------------------------------------------
!!> @brief Compute the area for each cell for affine mappings.
!!>
!!> If the mapping is affine, (e.g. Ax + b) then we don't need to sum
!!> the finer level areas.
!!>
!!> @param[in] mx, my the number of cells in the x and y directions
!!> @param[in] mbc the number of ghost cells
!!> @param[in] dx, dy the spacings in the x and y direcitons
!!> @param[in] xlower, ylower the lower left coordinate of the patch
!!> @param[in] blockno the block number
!!> @param[out] area the area of each cell
!!> @param[in] quadsize the length of the quad
!!> @param[in] quadstore stores a group of cell values
!!> @param[in] ghost_only
!! ---------------------------------------------------------
 
subroutine fclaw3d_metric_fort_compute_volume_affine(mx,my,mz, mbc,xd,yd,zd, & 
          volume, faceareas, ghost_only)
    implicit none

    integer mx,my,mz, mbc
    integer ghost_only

    double precision volume(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    double precision faceareas(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3)

    double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)

    integer i,j, k, m, icell, jcell, kcell
    double precision hex(0:1,0:1,0:1,3)
    double precision hex_compute_volume
    logical hex_is_volume_interior

    double precision area(3)

    do j = -mbc,my+mbc+1
        do i = -mbc,mx+mbc+1
            do k = -mbc,mz+mbc+1
                if (hex_is_volume_interior(mx,my,mz,i,j,k) .and. & 
                    ghost_only .eq. 1) then
                    cycle
                endif
                do icell = 0,1
                    do jcell = 0,1
                        do kcell = 0,1
                            hex(icell,jcell,kcell,1) = xd(i+icell,j+jcell,k+kcell)
                            hex(icell,jcell,kcell,2) = yd(i+icell,j+jcell,k+kcell)
                            hex(icell,jcell,kcell,3) = zd(i+icell,j+jcell,k+kcell)
                        end do
                    end do
                end do

                volume(i,j,k) = hex_compute_volume(hex)

                call hex_compute_surf_area(hex,area)
                do m = 1,3
                    faceareas(i,j,k,m) = area(m)
                end do
            end do
        end do
    end do
end subroutine fclaw3d_metric_fort_compute_volume_affine

!! ---------------------------------------------------------
!!> @brief Check if the index is an interior index
!!>
!!> @param[in] mx, my, mz the number of cell sin the x and y directions
!!> @param[in] i, j, k the index
!!> @return true if the index is interior
!! ---------------------------------------------------------
logical function hex_is_volume_interior(mx,my,mz,i,j,k)
    implicit none
    integer mx,my,mz,i,j,k

    hex_is_volume_interior = i .ge. 0 .and. i .le. mx+1 .and. & 
             j .ge. 0 .and. j .le. my+1 .and. & 
             k .ge. 0 .and. k .le. mz+1 
end function hex_is_volume_interior

!! Since basis vectors are computed on a coarse grid (assuming coarse rule
!! surface approximation at coarse grid cell face), we can compute a basis
!! independent of the volume calculation.  This is needed, for example, if
!! we re-build a coarse grid from fine grids
subroutine fclaw3d_metric_fort_compute_basis(mx, my, mz, mbc, xd, yd, zd, & 
          xrot, yrot, zrot, ghost_only)
    implicit none

    integer mx,my,mz, mbc
    integer ghost_only

    double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)
    double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2,-mbc:mz+mbc+2)

    double precision xrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
    double precision yrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
    double precision zrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)

    integer i,j, k, ii,jj, icell, jcell, kcell
    !! double precision xcorner, ycorner, zcorner
    !! double precision xe, ye,ze, xp1, yp1, zp1
    double precision hex(0:1,0:1,0:1,3)

    logical hex_is_volume_interior

    double precision rot(3,3,3)

    integer*8 map_context_ptr, fclaw_map_get_context

    map_context_ptr = fclaw_map_get_context()

    do j = -mbc,my+mbc+1
        do i = -mbc,mx+mbc+1
            do k = -mbc,mz+mbc+1
                if (hex_is_volume_interior(mx,my,mz,i,j,k) .and. & 
                    ghost_only .eq. 1) then
                    cycle
                endif
                !! xe = xlower + (i-1)*dx
                !! ye = ylower + (j-1)*dy
                !! ze = zlower + (k-1)*dz
                do icell = 0,1
                    do jcell = 0,1
                        do kcell = 0,1
!!                            xcorner = xe + icell*dx
!!                            ycorner = ye + jcell*dy
!!                            zcorner = ze + kcell*dz
!!                            call fclaw3d_map_c2m(map_context_ptr, & 
!!                                blockno,xcorner,ycorner,zcorner,xp1,yp1,zp1)
                            hex(icell,jcell,kcell,1) = xd(i+icell,j+jcell,k+kcell)
                            hex(icell,jcell,kcell,2) = yd(i+icell,j+jcell,k+kcell)
                            hex(icell,jcell,kcell,3) = zd(i+icell,j+jcell,k+kcell)
                        end do
                    end do
                end do
                call hex_compute_basis(hex,rot)
                do jj = 1,3
                    do ii = 1,3
                        xrot(i,j,k,ii,jj) = rot(1,ii,jj)
                        yrot(i,j,k,ii,jj) = rot(2,ii,jj)
                        zrot(i,j,k,ii,jj) = rot(3,ii,jj)
                    end do
                end do
            end do
        end do
    end do
end subroutine fclaw3d_metric_fort_compute_basis


subroutine hex_compute_basis(hex,rot)
    implicit none

    double precision hex(0:1,0:1,0:1,3)

    double precision Jb(3,3), Jinv(3,3), xcv(3), ycv(3), zcv(3)
    double precision eye(3,3), rot(3,3,3)
    double precision v(3,3), alpha(3), rhs(3,3)
    integer IPIV(3)

    double precision z000(3),z100(3),z010(3),z001(3), z110(3), & 
           z101(3), z011(3), z111(3)

    double precision a000(3),a100(3),a010(3),a001(3), a110(3), & 
           a101(3), a011(3), a111(3)

    logical is_id
    integer i,j,k, n, ni, info
    double precision sum


    do i = 1,3
        do j = 1,3
            if (i == j) then
                eye(i,j) = 1.d0
            else
                eye(i,j) = 0.d0
            endif
        end do
    end do


    !! Get centers of faces so we can find basis vectors at each face.
    do i = 1,3
        xcv(i) = 0.5d0
        ycv(i) = 0.5d0
        zcv(i) = 0.5d0
    end do

    xcv(1) = 0.d0
    ycv(2) = 0.d0
    zcv(3) = 0.d0

    !! # Make notation easy..
    do i = 1,3
        z000(i) = hex(0,0,0,i)
        z100(i) = hex(1,0,0,i)
        z010(i) = hex(0,1,0,i)
        z001(i) = hex(0,0,1,i)
        z110(i) = hex(1,1,0,i)
        z101(i) = hex(1,0,1,i)
        z011(i) = hex(0,1,1,i)
        z111(i) = hex(1,1,1,i)
    end do

    !! # Get coefficients for trilinear map.
    do i = 1,3
        a000(i) = z000(i)
        a001(i) = z001(i) - z000(i)
        a010(i) = z010(i) - z000(i)
        a100(i) = z100(i) - z000(i)

        a111(i) = z111(i) - z110(i) - z101(i) - & 
              z011(i) + z100(i) + z010(i) + z001(i) - z000(i)
        a011(i) = z011(i) - z010(i) - z001(i) + z000(i)
        a101(i) = z101(i) - z100(i) - z001(i) + z000(i)
        a110(i) = z110(i) - z100(i) - z010(i) + z000(i)
    end do

    !! # Start computing basis vectors.
    do n = 1,3  !! Loop over each face.

        do i = 1,3 !! Get three rows of Jacobian
            Jb(i,1) = a100(i) + a110(i)*ycv(n)  + a101(i)*zcv(n) + & 
                 a111(i)*ycv(n)*zcv(n)
            Jb(i,2) = a010(i) + a110(i)*xcv(n)  + a011(i)*zcv(n) + & 
                 a111(i)*xcv(n)*zcv(n)
            Jb(i,3) = a001(i) + a101(i)*xcv(n)  + a011(i)*ycv(n) + & 
                 a111(i)*xcv(n)*ycv(n)
        end do

        do i = 1,3
            IPIV(i) = 0
            do j = 1,3
                rhs(i,j) = eye(i,j)
            end do
        end do

        !! # Compute inverse of Jacobian.
        !!  call dgesv(3,3,Jb,3,IPIV,rhs,3,info)
        call hex_compute_Jinv(Jb,rhs,info)

         if (info .ne. 0) then
            !! # The Jacobian is singular;  set rotation to identity and hope
            !! # that it is never really needed in other clawpack calculations...
            do i = 1,3
                !! For Clawpack routines, frame is 
                !! ni = mod(i+n-2,3) + 1
                ni = i
                do j = 1,3
                    rot(n,i,j) = eye(ni,j)
                end do
            end do

            !! # Skip to next face (cycle through n loop)
            cycle

        endif

        !! # Now we have to permute Jinv to get normal to face n in first row.
        !! # We then orthogonalize other rows against this vector.
        do i = 1,3
            ni = mod(i+n-2,3) + 1
            do j = 1,3
                Jinv(i,j) = rhs(ni,j)
            end do
            write(10,'(3F16.8)') (Jinv(i,j),j=1,3)
        end do
        write(10,*) ' '

        !! # Do Gram-Schmidt on rows of Jinv, to get orthogonal basis.
        do i = 1,3

            !! # Project onto each previous vector v(i,:)
            alpha(1) = 0
            do k = 1,i-1
                alpha(k) = 0
                do j = 1,3
                    alpha(k) = alpha(k) + Jinv(i,j)*v(k,j)
                end do
            end do

            sum = 0
            do j = 1,3
                v(i,j) = Jinv(i,j)
                do k = 1,i-1
                    v(i,j) = v(i,j) - alpha(k)*v(k,j)
                end do
                sum = sum + v(i,j)*v(i,j)
            end do

            do j = 1,3
                v(i,j) = v(i,j)/dsqrt(sum)
            end do
        end do


        !! # We now have our three basis vectors for face n. Now, permute them
        !! back so that row 1 is always normal in direction xi, row 2 in
        !! direction eta, and row3 in direction zeta.

        !! # If we don't do this, then u is always the direction normal to the
        !! current face, v is always the first transverse, w the second.  The
        !! Riemann solvers as they are currently written assume that this is the
        !! case.

        !! ## Test to see if we are storing the identity matrix;  if so, we
        !! shouldn't store it above.
        !! is_id = .true.
         do i = 1,3
            !! ni = mod(i+n-2,3) + 1
            ni = i
            do j = 1,3
                if (abs(v(i,j) - eye(i,j)) > 1e-12) then
                    is_id = .false.
                endif
                rot(n,ni,j) = v(i,j)
            end do
        end do
    end do

end subroutine hex_compute_basis


!! ---------------------------------------------------------
!!> @brief Approximate the volume based on the corners of the mesh cell
!!>
!!> This is the volume element based on a trilinear approximation to
!!> the hexahedral cell defined by the corners of the mesh cell.
!!>
!!> @param[in] hex the coordinates of the cell nodes
!!> @return the exact volume of the hexahedral cell
!! ---------------------------------------------------------



double precision function hex_compute_volume(hex)
    implicit none


    double precision hex_dot_cross, volume
    double precision hex(0:1,0:1,0:1,3)
    double precision u(3),v(3),w(3), p(3), q(3), r(3)

    integer i,j,k,ip,jp,kp,m, ii, jj, kk

    do j = 1,3
        u(j) = hex(0,1,1,j) - hex(0,0,0,j)
        v(j) = hex(1,0,1,j) - hex(0,0,0,j)
        w(j) = hex(1,1,0,j) - hex(0,0,0,j)
    end do
    volume = dabs(hex_dot_cross(u,v,w))/6.d0

    do j = 1,3
        u(j) = hex(0,0,1,j) - hex(1,1,1,j)
        v(j) = hex(0,1,0,j) - hex(1,1,1,j)
        w(j) = hex(1,0,0,j) - hex(1,1,1,j)
    end do
    volume = volume + dabs(hex_dot_cross(u,v,w))/6.d0

    do i = 0,1
        do j = 0,1
            do k = 0,1
                ip = mod(i+1,2)
                jp = mod(j+1,2)
                kp = mod(k+1,2)

                do m = 1,3
                    p(m) = hex(ip,j,k,m) - hex(i,j,k,m)
                    q(m) = hex(i,jp,k,m) - hex(i,j,k,m)
                    r(m) = hex(i,j,kp,m) - hex(i,j,k,m)
                end do

                volume = volume + dabs(hex_dot_cross(p,q,r))/6.d0
            end do
        end do
    end do
    volume = volume/2.d0

    if (volume .eq. 0) then
        do ii = 0,1
            do jj = 0,1
                do kk = 0,1
                    write(6,*) hex(ii,jj,kk,1)
                    write(6,*) hex(ii,jj,kk,2)
                    write(6,*) hex(ii,jj,kk,3)
                end do
            end do
        end do
        write(6,*) 'volume is zero'
        stop
    endif

    hex_compute_volume = volume
end function hex_compute_volume

double precision function hex_dot_cross(u,v,w)
    implicit none
    double precision u(3), v(3), w(3)

    hex_dot_cross = u(1)*(v(2)*w(3) - v(3)*w(2)) - & 
                u(2)*(v(1)*w(3) - v(3)*w(1)) + & 
                u(3)*(v(1)*w(2) - v(2)*w(1))

    return
end function hex_dot_cross


!! # This will give us the exact surface area of each face in the case
!! # the surface is flat.  Otherwise, this seems like it should give a
!! # reasonable answer, but this is merely a conjecture!

subroutine hex_compute_surf_area(hex,area)
    implicit none

    double precision hex(0:1,0:1,0:1,3)
    double precision w(3), p(3), q(3), area(3)
    double precision quad(0:1,0:1,3)

    integer n,i,j,ip,jp,m
    double precision d

    do n = 1,3
        do i = 0,1
            do j = 0,1
                do m = 1,3
                    if (n == 1) then
                        quad(i,j,m) = hex(0,i,j,m)
                    elseif (n == 2) then
                        quad(i,j,m) = hex(i,0,j,m)
                    else
                        quad(i,j,m) = hex(i,j,0,m)
                    endif
                end do
            end do
        end do

        !! area of face is length of cross product of vectors on 
        !! ruled surface.  
        area(n) = 0
        do i = 0,1
            do j = 0,1
                ip = mod(i+1,2)
                jp = mod(j+1,2)

                do m = 1,3
                    p(m) = quad(ip,j,m) - quad(i,j,m)
                    q(m) = quad(i,jp,m) - quad(i,j,m)
                enddo
                w(1) =   p(2)*q(3) - p(3)*q(2)
                w(2) = -(p(1)*q(3) - p(3)*q(1))
                w(3) =   p(1)*q(2) - p(2)*q(1)

                d = w(1)*w(1) + w(2)*w(2) + w(3)*w(3)
                area(n) = area(n) + dsqrt(d)/2.d0
            end do
        end do
        area(n) = area(n)/2.d0
    enddo !! face loop
end subroutine hex_compute_surf_area


!! # This computes the inverse of a 3x3 matrix using Cramer's rule
!! # Note : Only if det(Jac) == 0 (exactly) will it report that
!! # the matrix is singular (i.e. info == 1).   It will not detect
!! # potential ill-conditioning.
subroutine hex_compute_Jinv(Jac,Jinv,info)
    implicit none

    double precision Jac(3,3), Jinv(3,3)
    double precision hex_dot_cross, detJ, s
    integer i,j, k(3,2), info

    data k /2, 1, 1, 3, 3, 2/

    info = 0

    !! # Compute determinant of Jacobian
    detJ  = hex_dot_cross(jac(1,1),jac(1,2),jac(1,3))
    if (detJ .eq. 0) then
        info = 1
        return
    endif

    !! # Apply Cramer's rule to get inverse.
    s = 1.d0/detJ
    do j = 1,3
        do i = 1,3
            Jinv(i,j) = s*(Jac(k(j,1),k(i,1))*Jac(k(j,2),k(i,2)) - & 
                           Jac(k(j,1),k(i,2))*Jac(k(j,2),k(i,1)))
            s = -s
        end do
    end do
end subroutine hex_compute_Jinv


!! ---------------------------------------------------------
!!> @brief Averages the volume from a fine grid to a coarse grid (3dx only). 
!!>        
!!>
!!> @param[in] mx, my the number of cells in the x and y directions
!!> @param[in] mbc the number of ghost cells
!!> @param[in] volfine the cell areas of the fine patch
!!> @param[out] volcoarse the cell areas of the coarse patch
!!> @param[in] igrid the index of the fine grid in the child array
!! ---------------------------------------------------------
subroutine fclaw3dx_metric_fort_average_volume(mx,my,mz, mbc, & 
           volcoarse, volfine, igrid)
    implicit none

    integer mx,my,mz, mbc, igrid

    !! # these will be empty if we are not on a manifold.
    double precision volcoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
    double precision   volfine(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

    integer i,j, k, ig, jg, ic_add, jc_add, ii, jj
    double precision sum

    !! # This should be refratio*refratio.
    integer i1,j1, r2, m
    integer rr2
    parameter(rr2 = 4)
    integer i2(0:rr2-1),j2(0:rr2-1)
    double precision kf

    integer p4est_refineFactor, refratio


    p4est_refineFactor = 2
    refratio = 2

    !! # 'iface' is relative to the coarse grid

    r2 = refratio*refratio
    if (r2 .ne. rr2) then
        write(6,*) 'average_face_ghost (claw2d_utils.f): Refratio**2 is not equal to rr2'
        stop
    endif


    !! # Get (ig,jg) for grid from linear (igrid) coordinates
    ig = mod(igrid,refratio)
    jg = (igrid-ig)/refratio

    !! # Get rectangle in coarse grid for fine grid.
    ic_add = ig*mx/p4est_refineFactor
    jc_add = jg*my/p4est_refineFactor

    r2 = refratio*refratio
    !! This only works for extruded mesh version.  
    !! We include one layer of ghost cells on the coarse grid 
    do k = 1,mz
        do j = 0,my/p4est_refineFactor+1
            do i = 0,mx/p4est_refineFactor +1
                i1 = i+ic_add
                j1 = j+jc_add
                m = 0
                do jj = 1,refratio
                    do ii = 1,refratio
                        i2(m) = (i-1)*refratio + ii
                        j2(m) = (j-1)*refratio + jj
                        m = m + 1
                    end do
                end do
                sum = 0
                do m = 0,r2-1
                    kf = volfine(i2(m),j2(m),k)
                    sum = sum + kf
                end do
                volcoarse(i1,j1,k) = sum
            end do
        end do
    end do

    !! # Compute area in the ghost cells

end subroutine fclaw3dx_metric_fort_average_volume


!! ---------------------------------------------------------
!!> @brief Averages the volume from a fine grid to a coarse grid (3dx only). 
!!>        
!!>
!!> @param[in] mx, my the number of cells in the x and y directions
!!> @param[in] mbc the number of ghost cells
!!> @param[in] volfine the cell areas of the fine patch
!!> @param[out] volcoarse the cell areas of the coarse patch
!!> @param[in] igrid the index of the fine grid in the child array
!! ---------------------------------------------------------
subroutine fclaw3dx_metric_fort_average_facearea(mx,my,mz, mbc, & 
           fa_coarse, fa_fine, igrid)
    implicit none

    integer mx,my,mz, mbc, igrid

    !! # these will be empty if we are not on a manifold.
    double precision fa_coarse(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3)
    double precision   fa_fine(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3)

    integer i,j, k, ig, jg, ic_add, jc_add, ii, jj
    double precision sum

    !! # This should be refratio*refratio.
    integer i1,j1, r2, m
    integer rr2
    parameter(rr2 = 4)
    integer i2(0:rr2-1),j2(0:rr2-1)
    double precision kf

    integer p4est_refineFactor, refratio


    p4est_refineFactor = 2
    refratio = 2

    !! # 'iface' is relative to the coarse grid

    r2 = refratio*refratio
    if (r2 .ne. rr2) then
        write(6,*) 'average_face_ghost (claw2d_utils.f): Refratio**2 is not equal to rr2'
        stop
    endif


    !! # Get (ig,jg) for grid from linear (igrid) coordinates
    ig = mod(igrid,refratio)
    jg = (igrid-ig)/refratio

    !! # Get rectangle in coarse grid for fine grid.
    ic_add = ig*mx/p4est_refineFactor
    jc_add = jg*my/p4est_refineFactor

    r2 = refratio*refratio
    !! This only works for extruded mesh version.
    do k = 1,mz
        do j = 0,my/2+1
            do i = 0,mx/2+1
                i1 = i+ic_add
                j1 = j+jc_add
                m = 0

                !! Set up linear indexing for 2x2 quad grid arrangement
                do jj = 1,2
                    do ii = 1,2
                        i2(m) = 2*(i-1) + ii
                        j2(m) = 2*(j-1) + jj
                        m = m + 1
                    end do
                end do

                !! Area of x-face
                if (igrid .eq. 0 .or. igrid .eq. 2) then
                    fa_coarse(i1,j1,k,1) = fa_fine(i2(0),j2(0),k,1) + & 
                                           fa_fine(i2(2),j2(2),k,1)
                endif

                !! Area of y-face
                if (igrid .eq. 0 .or. igrid .eq. 1) then
                    fa_coarse(i1,j1,k,2) = fa_fine(i2(0),j2(0),k,2) + & 
                                           fa_fine(i2(1),j2(1),k,2)
                endif

                !! Area of z-face
                sum = 0
                do m = 0,3
                    kf = fa_fine(i2(m),j2(m),k,3)
                    sum = sum + kf
                end do
                fa_coarse(i1,j1,k,3) = sum
            end do
        end do
    end do

    !! # Compute area in the ghost cells

end subroutine fclaw3dx_metric_fort_average_facearea
