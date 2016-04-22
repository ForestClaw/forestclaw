SUBROUTINE geoclaw_setup()

  USE geoclaw_module, ONLY: set_geo
  USE qinit_module, ONLY: set_qinit

  IMPLICIT NONE

  CALL set_geo()                    !# sets basic parameters g and coord system
  CALL set_qinit()                  !# specifies file with dh if this used instead

  END SUBROUTINE geoclaw_setup
