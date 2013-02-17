SUBROUTINE set_modules_defaults()
  USE geoclaw_module
  USE topo_module
  implicit none

  !! These were being set in amr2ez_geo.f, so thought I should set things to
  !! 0 here.  This is called from amr_geoclaw_setup_problem().

  mgauges = 0
  mtopofiles = 0
  mregions = 0
  idtopo = 0
  iqinit = 0
  icoordsys = 0
  coeffmanning = 0.d0
  frictiondepth = 0.d0

END SUBROUTINE set_defaults
