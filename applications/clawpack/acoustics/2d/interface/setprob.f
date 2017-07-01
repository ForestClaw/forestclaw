      subroutine interface_setprob(rhol,cl,rhor,cr)
      implicit none
      double precision rhol,cl,rhor,cr
      double precision rhol_com,cl_com,rhor_com,cr_com
      common /comaux/ rhol_com,cl_com,rhor_com,cr_com

c     # this should be set from user options.
      rhol_com = rhol
      cl_com = cl
      rhor_com = rhor
      cr_com = cr

      return
      end
