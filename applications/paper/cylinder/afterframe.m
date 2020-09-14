setviews;

[~,~,~,R,H,~] = read_vars();

fprintf('%6s %16.8e\n','qmin',qmin);
fprintf('%6s %16.8e\n\n','qmax',qmax);

conservation_cylinder = false;
if (conservation_cylinder)
    th = 2*pi*(0.25 + 1.d0/32.d0);
    x0 = cos(th);
    y0 = sin(th);
    z0 = alpha + 0.01;
    hold on;
    plot3(x0,y0,z0,'k.','markersize',30);
    plot3(-x0,y0,z0,'k.','markersize',30);
    hold off;
end

% yrbcolormap;
colormap(parula);

caxis([-0.2,1]);
if (qmin < qmax)
    caxis([qmin,qmax]);
end

showpatchborders;
setpatchborderprops('linewidth',1);

hc = colorbar;
set(hc,'fontsize',16);

daspect([1,1,1]);

axis off;

set(gcf,'clipping','off')
set(gca,'clipping','off');

prt = false;
NoQuery = false;
if (prt)
    delete(get(gca,'title'));
    if (mq == 1)
        fname = sprintf('soln_torus%02d.png',Frame);
    else
        hidepatchborders(6);
        fname = sprintf('error_torus%02d.png',Frame);
    end
    print('-dpng','-r240',fname);
end
shg

clear afterframe;
clear mapc2m;
