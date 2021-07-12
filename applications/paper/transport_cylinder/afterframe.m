setviews;

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

showgridlines;

colormap(parula);
showpatchborders;
setpatchborderprops('linewidth',1);
daspect([1,1,1]);
colorbar;

caxis([-0.2,1])
if (abs(qmin-qmax) < 1e-8)
    if (qmin < qmax)
        caxis([qmin,qmax]);
    end
end
% caxis([-1,1]*1e-8);
axis off;

view(3)
set(gca,'clipping','off');

prt = false;
NoQuery = false;
if (prt)
    delete(get(gca,'title'));
    fname = sprintf('torus%02d.png',Frame);
    print('-dpng','-r640',fname);
end
shg
