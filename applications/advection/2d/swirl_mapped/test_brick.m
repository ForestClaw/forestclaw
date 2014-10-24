function test_brick(blockno)

close all;

set_blocknumber(blockno);

mx = 8;
mbc = 0;

minlevel = 2;
plot_coarse = true;


% Plot coarse grid
if (plot_coarse)
    mxc = mx*2^minlevel;
    mbcc = mbc;
    dxc = 1/mxc;
    x = linspace(-mbcc*dxc,1+mbcc*dxc,mxc+2*mbcc+1);
    y = x;
    [xm,ym] = meshgrid(x,y);

    [xp,yp,zp] = mapc2m(xm,ym);
        
    ph = patch(surf2patch(xp,yp,zp));
    set(ph,'edgecolor','r');
    set(ph,'facecolor','none');
end

daspect([1 1 1])


shg

end


        