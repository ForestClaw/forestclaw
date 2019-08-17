setviews;

plot_contour = false;
if plot_contour
    yrbcolormap;
    hold on;
    % assume we are doing the Gaussian problem with error
    [xp,yp,zp] = torus_soln_contour(t,alpha,R,period);
    plot3(xp,yp,zp,'k','linewidth',2);
    hold off;
    hidepatchborders;
    view(3)
    % camlight;
    axis off
end

fprintf('%6s %16.8e\n','qmin',qmin);
fprintf('%6s %16.8e\n\n','qmax',qmax);

% yrbcolormap;
colormap(parula);

showpatchborders;
setpatchborderprops('linewidth',1);

daspect([1,1,1]);

caxis([0, 1]);
colorbar
axis off;

view(3);
view(vtop);

if (mq >= 3)
    caxis([qmin,qmax]);
    colorbar;
end

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
