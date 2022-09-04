      subroutine fclaw_map_set_context(context)
      implicit none

      integer*8 context, context_com
      common /comcontext/ context_com

      context_com = context
      end

      integer*8 function fclaw_map_get_context()
      implicit none

      integer*8 context_com
      common /comcontext/ context_com

      fclaw_map_get_context = context_com
      return
      end
