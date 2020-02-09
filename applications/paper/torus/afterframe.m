setviews;

[~,~,~,alpha,beta,~] = read_vars();

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

conservation_torus = false;
if (conservation_torus)
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

if (mq == 3)
    caxis([-1,1]*1e-3);
else
    if Frame == 0
        caxis([-0.3,1.0]);
    else
        if (mq == 1)
            caxis([-0.7,2.1])    
        else
            caxis([-1,1]*1e-3);
        end
    end
end

showpatchborders;
setpatchborderprops('linewidth',1);

hc = colorbar;
set(hc,'fontsize',16);

daspect([1,1,1]);

axis off;

set(gcf,'clipping','off')
set(gca,'clipping','off');

view([0,27]);

prt = true;
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
