yrbcolormap
axis([0 1 0 1 0 1])
daspect([1 1 1]);


fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);


showpatchborders;
setpatchborderprops('linewidth',1);
setpatchbordercolor('k');

showcubes;
showslices;

caxis([0,1]);

h = surflight;

axis off;

shg;

clear afterframe;
