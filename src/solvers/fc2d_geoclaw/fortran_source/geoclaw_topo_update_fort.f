      subroutine fc2d_geoclaw_topo_update(t)
      use topo_module, only: topo_finalized

      if (.not. topo_finalized) then
         call topo_update(t)
      endif
      end