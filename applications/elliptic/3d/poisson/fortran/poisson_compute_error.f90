subroutine poisson_compute_error(blockno, mx,my,mz,mbc,mfields, &
   dx,dy,dz,xlower,ylower,zlower,t,rhs,error,soln)
  implicit none

  integer :: mx,my,mz,mbc,mfields, blockno
  double precision :: dx, dy, dz, xlower, ylower, zlower, t
  double precision, dimension(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,mfields) :: rhs
  double precision, dimension(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,mfields) :: error
  double precision, dimension(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,mfields) :: soln

  integer :: i,j,k,m
  double precision :: xc,yc,zc, poisson_qexact

  ! Assume a single field variable only
  do k = 1,mz
      do j = 1,my
           do i = 1,mx
               zc = zlower + (k-0.5)*dz
               yc = ylower + (j-0.5)*dy
               xc = xlower + (i-0.5)*dx

               soln(i,j,k,1) = poisson_qexact(xc,yc,zc)
               do m = 1,mfields
                    error(i,j,k,m) = rhs(i,j,k,m) - soln(i,j,k,1)
               end do
           end do
      end do
  end do

end subroutine poisson_compute_error
