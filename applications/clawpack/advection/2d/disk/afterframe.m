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
plot(xout,yout,'k','linewidth',2);
% fprintf('Area of filament %24.16f\n',polyarea(xout,yout));
hold off;

showpatchborders;
setpatchborderprops('linewidth',1)

view(2);
axis off
shg;


NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 41;
  filename = sprintf('filament_%04d.png',Frame)
  print('-dpng',filename);
end
