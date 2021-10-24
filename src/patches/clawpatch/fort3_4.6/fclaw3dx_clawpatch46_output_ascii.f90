subroutine fclaw3dx_clawpatch46_fort_header_ascii (matname1, &
    matname2, time,meqn,maux,ngrids)
    implicit none

    integer meqn,ngrids, maux

    character*11 matname1
    character*11 matname2
    double precision time
    integer matunit1, matunit2

    matunit1 = 10
    matunit2 = 15

    open(unit=matunit2,file=matname2)
    write(matunit2,1000) time,meqn,ngrids,maux,3
 1000 format(e30.20,'    time', /,  & 
           i5,'                 meqn'/,  & 
           i5,'                 ngrids'/, & 
           i5,'                 num_aux'/, & 
           i5,'                 num_dim')

    close(matunit2)

    open(unit=matunit1,file=matname1,status='replace')
    close(matunit1)

end subroutine  fclaw3dx_clawpatch46_fort_header_ascii

subroutine fclaw3dx_clawpatch46_fort_output_ascii(matname1, & 
    mx,my,mz,meqn,mbc, xlower,ylower, zlower, dx,dy,dz, & 
    q,patch_num,level,blockno,mpirank)

    implicit none

    character(len=11) matname1
    integer meqn,mbc,mx,my, mz
    integer patch_num
    integer level, blockno, mpirank
    double precision xlower, ylower,zlower, dx,dy, dz

    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer matunit1
    integer i,j,k,mq

    matunit1 = 10
    open(matunit1,file=matname1,position='append');

    call fclaw3dx_clawpatch46_fort_write_grid_header(matunit1, & 
           mx,my,mz,xlower,ylower, zlower, dx,dy,dz, patch_num,level, & 
           blockno,mpirank)


    if (meqn .gt. 5) then
        write(6,'(A,A,A,I5,A)')     & 
              'Warning (fclaw2d_fort3_write_grid_header.f) ', & 
              ': meqn > 5; change format statement 120.', & 
              '(meqn = ',meqn,')'
        stop
    endif

    do k = 1,mz
        do j = 1,my
            do i = 1,mx
                do mq = 1,meqn
                    if (abs(q(i,j,k,mq)) .lt. 1d-99) then
                        q(i,j,k,mq) = 0.d0
                    endif
                enddo
                write(matunit1,120) (q(i,j,k,mq),mq=1,meqn)
            end do
            write(matunit1,*) ' '
        end do
        write(matunit1,*) ' '
        !! # This statement is checked above (meqn <= 5)
    end do
  120 format (5E26.16)

    close(matunit1)

end subroutine fclaw3dx_clawpatch46_fort_output_ascii


subroutine fclaw3dx_clawpatch46_fort_write_grid_header (matunit1, & 
    mx,my,mz, xlower,ylower, zlower, dx,dy,dz, patch_num,level, & 
     blockno,mpirank)

    implicit none

    integer :: matunit1, mx, my, mz
    integer :: patch_num, level, blockno, mpirank
    double precision :: xlower, ylower,zlower, dx,dy, dz


    write(matunit1,1001) patch_num, level, blockno, mpirank, mx, my, mz
 1001 format(i5,'                 grid_number',/,  & 
             i5,'                 AMR_level',/,    & 
             i5,'                 block_number',/, & 
             i5,'                 mpi_rank',/,     & 
             i5,'                 mx',/,           & 
             i5,'                 my',/,           & 
             i5,'                 mz',/)


      write(matunit1,1002) xlower,ylower,zlower, dx,dy, dz
 1002 format(e24.16,'    xlow', /,  & 
            e24.16,'    ylow', /,   & 
            e24.16,'    zlow', /,   & 
            e24.16,'    dx', /,     & 
            e24.16,'    dy',/,      & 
            e24.16,'    dz',/)
end subroutine fclaw3dx_clawpatch46_fort_write_grid_header
