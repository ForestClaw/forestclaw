yrbcolormap
axis([0 2 0 2 0 1])
daspect([1 1 1]);


fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);


showpatchborders;
setpatchborderprops('linewidth',1);
setpatchbordercolor('k');

showcubes(5);
showslices;
hidesurfs
view(2);


caxis([0,1]);

h = surflight;

set(gca,'box','on');

%axis off;
axis on;

shg;

clear afterframe;
