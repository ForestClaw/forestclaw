subroutine get_aux_locations_n(ixyz,mcapa,locrot,locarea)
    implicit none

    integer ixyz,mcapa,locrot,locarea

    locarea = mcapa + ixyz
    locrot = mcapa + 4 + (ixyz-1)*9

    return
end subroutine get_aux_locations_n


!! # This gets the aux locations needed for the transverse
!! # Riemann solves
subroutine get_aux_locations_t(ixyz,icoor,mcapa,locrot,locarea)
    implicit none

    integer ixyz, icoor, locrot, locarea, mcapa
    integer nxyz

    nxyz = 0
    if (ixyz .eq. 1) then
        if (icoor .eq. 2) then
            !! # transverse solve is y-like
            !! # y-like direction is a y-face
            nxyz = 2
        elseif (icoor .eq. 3) then
            !! # transverse solve is z-like
            !! # z-like direction is a z-face
            nxyz = 3
        endif
    elseif (ixyz .eq. 2) then
        if (icoor .eq. 2) then
            !! # transverse solve is y-like
            !! # y-like direction is a z-face
            nxyz = 3
        else
            !! # transverse solve is z-like
            !! # y-like direction is an x-face
            nxyz = 1
        endif
    elseif (ixyz .eq. 3) then
        if (icoor .eq. 2) then
            !! # transverse solve is y-like
            !! # y-like direction is an x-face
            nxyz = 1
        elseif (icoor .eq. 3) then
            !! # transverse solve is z-like
            !! # z-like direction is a y-face
            nxyz = 2
        endif
    endif

    locarea = mcapa + nxyz
    locrot = mcapa + 4 + (nxyz-1)*9

end subroutine get_aux_locations_t


!! # This gets the aux locations needed for the double transverse
!! # Riemann solves
subroutine get_aux_locations_tt(ixyz,icoor,mcapa,locrot,locarea,irot)
    implicit none

    integer ixyz, icoor, locrot, locarea,mcapa, irot
    integer nxyz

    nxyz = 0
    if (ixyz .eq. 1) then
        !! # (x-like, y-like, z-like) = (x,y,z)
        if (icoor .eq. 2) then
            !! # transverse solve is y-like solve
            nxyz = 2
        elseif (icoor .eq. 3) then
            !! # transverse solve is z-like solve
            nxyz = 3
        endif
    elseif (ixyz .eq. 2) then
        !! # (x-like, y-like, z-like) = (y,z,x)
        if (icoor .eq. 2) then
            !! # transverse solve is z-like solve
            nxyz = 3
        elseif (icoor .eq. 3) then
            !! # transverse solve is x-like solve
            nxyz = 1
        endif
    elseif (ixyz .eq. 3) then
        !! # (x-like, y-like, z-like) = (z,x,y)
        if (icoor .eq. 2) then
            !! # transverse solve is x-like solve
            nxyz = 1
        elseif (icoor .eq. 3) then
            !! # transverse solve is y-like solve
            nxyz = 2
        endif
    endif

    locarea = mcapa + nxyz
    locrot = mcapa + 4 + (nxyz-1)*9
    irot = nxyz

end subroutine get_aux_locations_tt


subroutine rotate3(rot,velcomps)
    implicit none

    double precision  velcomps(3), rot(9)
    double precision v1, v2, v3

    v1 = velcomps(1)
    v2 = velcomps(2)
    v3 = velcomps(3)

    velcomps(1) = rot(1)*v1 + rot(2)*v2 + rot(3)*v3
    velcomps(2) = rot(4)*v1 + rot(5)*v2 + rot(6)*v3
    velcomps(3) = rot(7)*v1 + rot(8)*v2 + rot(9)*v3
end subroutine rotate3

subroutine rotate3_tr(rot,velcomps)
    implicit none

    double precision velcomps(3),rot(9)
    double precision v1, v2, v3

    v1 = velcomps(1)
    v2 = velcomps(2)
    v3 = velcomps(3)

    velcomps(1) = rot(1)*v1 + rot(4)*v2 + rot(7)*v3
    velcomps(2) = rot(2)*v1 + rot(5)*v2 + rot(8)*v3
    velcomps(3) = rot(3)*v1 + rot(6)*v2 + rot(9)*v3
end subroutine rotate3_tr


!! subroutine compute_binormal(rot)
!!     implicit none
!! 
!!     double precision rot(9)
!!     double precision w1,w2,w3,r
!! 
!!     w1 =  rot(2)*rot(6) - rot(3)*rot(5)
!!     w2 = -rot(1)*rot(6) + rot(3)*rot(4)
!!     w3 =  rot(1)*rot(5) - rot(2)*rot(4)
!!     r = dsqrt(w1*w1 + w2*w2 + w3*w3)
!! 
!!     rot(7) = w1/r
!!     rot(8) = w2/r
!!     rot(9) = w3/r
!! 
!! end subroutine compute_binormal
