      subroutine cudaclaw5_fort_update_q(meqn, mx, my, mbc, maux,
     &                                   dtdx, dtdy,qold,fp,fm,
     &                                   gp, gm, mcapa)
      implicit none

      integer meqn, mx, my, maux, mbc, mcapa
      double precision dtdx, dtdy

      double precision qold(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision fp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

      integer m, i, j

      do m = 1,meqn
         do i = 1,mx
            do j = 1,my
               if (mcapa .eq. 0) then
c                 # no capa array.  Standard flux differencing:
                  qold(m,i,j) = qold(m,i,j)
     &                  - dtdx * (fm(m,i+1,j) - fp(m,i,j))
     &                  - dtdy * (gm(m,i,j+1) - gp(m,i,j))
               else
c                 # with capa array.
                  qold(m,i,j) = qold(m,i,j)
     &                  -(dtdx*(fm(m,i+1,j) - fp(m,i,j))
     &                  + dtdy*(gm(m,i,j+1) - gp(m,i,j)))/aux(mcapa,i,j)
               endif
            enddo
         enddo
      enddo

      end


