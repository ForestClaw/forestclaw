yrbcolormap
axis([0 1 0 1 0 1])
daspect([1 1 1]);

showslices;
hidecubes;

setpatchborderprops('linewidth',1)

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);


showpatchborders;

caxis([0,1]);

axis off;

shg;

clear afterframe;
