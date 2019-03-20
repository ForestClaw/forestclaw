setviews;

alpha = 0.4;
s = 1e-2;    
alim = [-1-alpha,1+alpha];
axis([alim alim]);
daspect([1 1 1]);
view(vtop)

yrbcolormap;

showpatchborders(1:10);
setpatchborderprops('linewidth',1)
caxis([0,1])

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

if (mq == 3)
    % Plot the error
    ca = [-max([qmin,qmax]),max([qmin,qmax])];
else    
    % Plot the solution
    yrbcolormap;
    ca = [0,1];
end

% caxis(ca);
caxis([-1,1]*1e-8)



%
NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 8;
  axis([0 1 0 1]);
  filename = sprintf('annulus_%04d.png',Frame)
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
