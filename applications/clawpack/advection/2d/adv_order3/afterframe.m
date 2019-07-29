daspect([1,1,1]);

colormap(parula);
caxis([-1,1]);

showpatchborders;

fprintf('%15s %24.16e\n','qmin',qmin);
fprintf('%15s %24.16e\n','qmax',qmax);

shg

clear afterframe
