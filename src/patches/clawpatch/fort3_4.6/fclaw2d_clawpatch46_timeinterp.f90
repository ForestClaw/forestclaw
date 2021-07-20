subroutine fclaw2d_clawpatch46_fort3_timeinterp(mx,my,mz, mbc,meqn, & 
        psize, qcurr,qlast,qinterp,alpha,ierror)
    implicit none

    integer :: mx,my,mz,mbc,meqn,psize, ierror
    double precision :: alpha
    double precision ::   qcurr(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::   qlast(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: qinterp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer :: i,j,k,m,mint,count,count_final

    !! # Number of interior layers to compute.  Since we only
    !! # interpolate from time-interpolated levels, we only
    !! # need two layers.  If we were averaging, we'd need four.
    mint = mbc
    ierror = 0
    count = 1

    !! # Time interpolate only to the interior cells.

    do m = 1,meqn
        do k = 1,mz
            !! # Face 0
            do j = 1,my-mint
                do i = 1,mint
                    qinterp(i,j,k, m) = qlast(i,j,k,m) + & 
                        alpha*(qcurr(i,j,k,m)-qlast(i,j,k,m))
                    count = count + 1
                end do
            end do


            !! # Face 2
            do j = 1,mint
                do i = mint+1,mx
                    qinterp(i,j,k,m) = qlast(i,j,k,m) + & 
                        alpha*(qcurr(i,j,k,m)-qlast(i,j,k,m))
                    count = count + 1
                end do
            end do
 
            !! # Face 1
            do j = mint+1,my
                do i = mx-mint+1,mx
                    qinterp(i,j,k, m) = qlast(i,j,k,m) + & 
                         alpha*(qcurr(i,j,k,m)-qlast(i,j,k,m))
                    count = count + 1
                end do
            end do

            !! # Face 3
            do j = my-mint+1,my
                do i = 1,mx-mint
                    qinterp(i,j,k,m) = qlast(i,j,k,m) + &
                         alpha*(qcurr(i,j,k,m)-qlast(i,j,k,m))
                    count = count + 1
                end do
            end do
        end do
    end do

    count_final = count-1
    if (count_final .ne. psize) then
       ierror = 2
    endif


end
