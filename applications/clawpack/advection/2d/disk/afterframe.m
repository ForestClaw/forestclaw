yrbcolormap;
s = 0.05;
axis([-s 2+s -s 2+s])
daspect([1 1 1]);

showpatchborders(1:9);
caxis([0 1]);

N = 500;
if (t > 0)
    [xout,yout] = filament_soln(N,t);
else
    th = linspace(0,2*pi,N+1);
    xout = 0.25*cos(th) + 0.5;
    yout = 0.25*sin(th) + 1;
end
hold on;
plot(xout,yout,'k','linewidth',2);
fprintf('Area of filament %24.16f\n',polyarea(xout,yout));
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

clear afterframe
clear mapc2m
clear mapc2m_squareddisk
clear mapc2m_pillowdisk
clear parallelpartitions
