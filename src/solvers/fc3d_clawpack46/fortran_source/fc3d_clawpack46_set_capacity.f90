subroutine clawpack46_set_capacity(mx,my,mz,mbc,dx,dy,dz, & 
                                   volume,mcapa,maux,aux)
    implicit none

    integer :: mbc, mx, my, mz, maux, mcapa
    double precision :: dx, dy, dz
    double precision :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc, maux)
    double precision :: volume(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

    integer :: i,j,k
    double precision :: dxdydz

    dxdydz = dx*dy*dz

    do k = 1-mbc,mz+mbc
        do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
                aux(i,j,k,mcapa) = volume(i,j,k)/dxdydz
            end do
        end do
    end do
end subroutine clawpack46_set_capacity
