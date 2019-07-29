!! This would generally be supplied by the user

SUBROUTINE  rpn2qad_flux(meqn,maux,idir,iface,q, auxvec,flux)

  IMPLICIT NONE

  INTEGER meqn,maux,idir, m, iface
  DOUBLE PRECISION q(meqn), flux(meqn)
  DOUBLE PRECISION auxvec(maux)

  INTEGER color_equation
  COMMON /eqn_comm/ color_equation

  INTEGER qad_mode
  COMMON /qad_comm/ qad_mode

  DOUBLE PRECISION urot, g
  INTEGER k

  IF (color_equation .ne. 0) THEN
      DO m = 1,meqn
          flux(m) = 0
      ENDDO
  ELSE
      !! Cell-centered velocity has been projected to normals on each face.
      urot = auxvec(2+iface)  

!!    Use this version if only 2 lengths (left, bottom) are stored in 
!!    each cell and aux variables have been identified as "xleft", "yleft", etc.      
      IF (qad_mode .le. 1) THEN
          !! This won't work, but illustrates behavior for original qad
          g = auxvec(6+2*idir) !! Edge length for iface = 0,1,2,3
      ELSE
          !! Use this version if all four edge lengths are available in each cell. 
          g = auxvec(6+iface)  !! Edge length for iface = 0,1,2,3
      ENDIF

      DO m = 1,meqn
         !! # Scaling done here (unlike in ForestClaw)    
         flux(m) = g*urot*q(m)
      ENDDO
  ENDIF

END SUBROUTINE rpn2qad_flux
