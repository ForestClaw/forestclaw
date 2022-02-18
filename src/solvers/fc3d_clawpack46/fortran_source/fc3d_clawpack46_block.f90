subroutine fc3d_clawpack46_set_block(blockno)
    implicit none

    integer :: blockno, blockno_com
    common /comblock/ blockno_com

    blockno_com = blockno
    end subroutine fc3d_clawpack46_set_block

integer function fc3d_clawpack46_get_block()
    implicit none

    integer :: blockno_com
    common /comblock/ blockno_com

    fc3d_clawpack46_get_block = blockno_com
    return
end function fc3d_clawpack46_get_block

subroutine fc3d_clawpack46_unset_block()
    implicit none

    integer :: blockno_com
    common /comblock/ blockno_com

    blockno_com = -1
end subroutine fc3d_clawpack46_unset_block
