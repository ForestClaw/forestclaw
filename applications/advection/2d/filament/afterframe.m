s = 0.05;
axis([-s 2+s -s 2+s])
daspect([1 1 1]);
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
hold off;

showpatchborders(1:7);
setpatchborderprops('linewidth',1)

view(2);
axis off
shg;


NoQuery = 1;
prt = true;
if (prt)
  MaxFrames = 32;
  filename = sprintf('filament%03d.png',Frame);
  fprintf('Printing %s ...\n',filename);
  print('-dpng',filename);
end


clear afterframe
clear mapc2m
