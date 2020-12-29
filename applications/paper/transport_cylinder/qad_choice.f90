SUBROUTINE qad(valbig,mitot,mjtot,nvar,  &
    svdflx,qc1d,lenbc,lratiox,lratioy,hx,hy,  &
    maux,auxbig,auxc1d,delt,mptr)

!!  USE amr_module, only : timemult, nghost, max1d, maxaux, mwaves,nestlevel,rnode, &
!!       auxtype, node, method
  IMPLICIT NONE

  INTEGER mitot, mjtot, nvar, lenbc, lratiox, lratioy, maux
  INTEGER mptr
  DOUBLE PRECISION hx,hy,delt

  DOUBLE PRECISION valbig(nvar,mitot,mjtot)
  DOUBLE PRECISION auxbig(maux,mitot,mjtot)
  DOUBLE PRECISION qc1d(nvar,lenbc)
  DOUBLE PRECISION svdflx(nvar,lenbc)
  DOUBLE PRECISION auxc1d(maux,lenbc)
      

  INTEGER qad_mode
  COMMON /qad_comm/ qad_mode

  IF (qad_mode .eq. 0) THEN
!!    # Don't include conservation fix  
      return
  ELSEIF (qad_mode .eq. 1) THEN
!!    # Original qad
      CALL qad_original(valbig,mitot,mjtot,nvar,  &
                     svdflx,qc1d,lenbc,lratiox,lratioy,hx,hy,  &
                     maux,auxbig,auxc1d,delt,mptr)
  ELSEIF (qad_mode .eq. 2) THEN         
      CALL qad_modified(valbig,mitot,mjtot,nvar,  &
                        svdflx,qc1d,lenbc,lratiox,lratioy,hx,hy,  &
                        maux,auxbig,auxc1d,delt,mptr)
  ELSEIF (qad_mode .eq. 3) then
!!    # New qad
      CALL qad_new(valbig,mitot,mjtot,nvar, &
                   svdflx,qc1d,lenbc,lratiox,lratioy,hx,hy, &
                   maux,auxbig,auxc1d,delt,mptr)   
  ENDIF

END SUBROUTINE qad    
         
