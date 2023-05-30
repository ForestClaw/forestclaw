subroutine euler3d_setaux_manifold(mbc,mx,my,mz, mcapa, & 
               xlower,ylower,zlower,dx,dy,dz,maux, & 
               aux,blockno, xrot, yrot, zrot, volume, & 
               faceareas)
    implicit none

    integer mbc, mx, my, mz, maux, mcapa
    integer blockno
    double precision dx,dy, dz, xlower, ylower, zlower
    double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc, 1-mbc:mz+mbc,maux)

    double precision xrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
    double precision yrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
    double precision zrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)

    double precision volume(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc)
    double precision faceareas(-mbc:mx+mbc+1,-mbc:my+mbc+1, -mbc:mz+mbc+2,3)

    integer i,j, k, ii,jj,count
    double precision dxdydz, dxdy, dxdz,dydz

    dxdydz = dx*dy*dz
    dxdy = dx*dy
    dxdz = dx*dz
    dydz = dy*dz

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            do k = 1-mbc,mz+mbc
                aux(i,j,k,mcapa) = volume(i,j,k)/dxdydz
                !!write(6,*) aux(i,j,k,mcapa)

                aux(i,j,k,mcapa+1) = faceareas(i,j,k,1)/dydz
                aux(i,j,k,mcapa+2) = faceareas(i,j,k,2)/dxdz
                aux(i,j,k,mcapa+3) = faceareas(i,j,k,3)/dxdy
                !!write(6,*) aux(i,j,k,mcapa+1)
                !!write(6,*) aux(i,j,k,mcapa+2)
                !!write(6,*) aux(i,j,k,mcapa+3)
                !!write(6,*) ' '

                count = 1
                do ii = 1,3
                    do jj = 1,3
                        !! Store by row
                        aux(i,j,k,mcapa+3+ 0+count) = xrot(i,j,k,ii,jj)
                        aux(i,j,k,mcapa+3+ 9+count) = yrot(i,j,k,ii,jj)
                        aux(i,j,k,mcapa+3+18+count) = zrot(i,j,k,ii,jj)
                        count = count + 1
                    end do
                end do
            end do
        end do
    end do
    !!stop

end subroutine euler3d_setaux_manifold

