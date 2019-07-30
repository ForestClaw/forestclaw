setviews;

fprintf('%6s %16.8e\n','qmin',qmin);
fprintf('%6s %16.8e\n\n','qmax',qmax);

% yrbcolormap;
colormap(parula);

% showgridlines(1:5);

showpatchborders;
% hidepatchborders;
setpatchborderprops('linewidth',1);

daspect([1,1,1]);

caxis([0, 1]);
axis off;
set(gca,'clipping','off');

view(3);

prt = false;
NoQuery = false;
if (prt)
    delete(get(gca,'title'));
    fname = sprintf('torus%02d.png',Frame);
    print('-dpng','-r640',fname);
end
shg

clear afterframe;
clear mapc2m;
