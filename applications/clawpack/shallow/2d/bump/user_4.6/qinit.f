      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Shallow water with radial dam break problem, h = hin inside
c     # circle specified in fdisc.f
c
       implicit double precision (a-h,o-z)
       dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

       integer blockno, fc2d_clawpack46_get_block

       common /comic/ hin,hout

       blockno = fc2d_clawpack46_get_block()

       do 20 i=1-mbc,mx+mbc
          xlow = xlower + (i-1.d0)*dx
          xc = xlower + (i-0.5d0)*dx
          do 20 j=1-mbc,my+mbc
             ylow = ylower + (j-1.d0)*dy
             yc = ylower + (j-0.5d0)*dy
             q(i,j,1) = 0.1d0 + exp(-200.d0*(xc**2 + yc**2))

c             call cellave2(blockno,xlow,ylow,dx,dy,win)
c             q(i,j,1) = hin*win + hout*(1.d0-win)
             q(i,j,2) = 0.d0
             q(i,j,3) = 0.d0
  20         continue
       return
       end
