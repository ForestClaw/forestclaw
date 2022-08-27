yrbcolormap;
s = 0.05;
ax = 0;
bx = 2;
ay = 0;
by = 2;
axis([ax-s bx+s ay-s by+s])
daspect([1 1 1]);

showpatchborders(1:9);
caxis([0 1]);

N = 500;
plot_soln = true;
if plot_soln
    if (t > 0)
        [xout,yout] = filament_soln(N,t);
    else
        R_init = 0.25;
        xc_init = 0.5;
        yc_init = 1.0;
        th = linspace(0,2*pi,N+1);
        xout = R_init*cos(th) + xc_init;
        yout = R_init*sin(th) + yc_init;
    end
    hold on;
    for m = 1:length(zSliceCoords)
        zout = zeros(size(xout)) + zSliceCoords(m);
        plot3(xout,yout,zout,'k','linewidth',3);
    end
    hold off;
end

showpatchborders;
setpatchborderprops('linewidth',1)
showsurfs;

axis([0 2 0 2 0 1])
daspect([1 1 0.25])
view(2);
hidesurfs;
% showcubes(4);
axis off



NoQuery = 0;
prt = false;
if (prt)
    MaxFrames = 41;
    filename = sprintf('filament_%04d.png',Frame)
    print('-dpng',filename);
end

shg
