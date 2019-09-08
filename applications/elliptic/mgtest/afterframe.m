ax = -1;
bx = 1;
ay = -1;
by = 1;
s = 1e-2;
axis([ax-s bx+s ay-s by+s])
daspect([1 1 1]);
% axis off;

showpatchborders;
setpatchborderprops('linewidth',1);

view(2);

if (length(amrdata) == 1)
    % We are only on a single grid
    figure(2);
    clf;    
    h = surf(xcenter,ycenter,q);
    set(h,'edgecolor','none');
    view(3);
    axis square;
    set(gcf,'color','k');
    set(gcf,'clipping','off');
    axis off;
    camlight;
    figure(1);
end
if Frame == 0 && mq == 1
    caxis([-1,1]*1e-3);
end
if Frame == 1 && mq == 1
    caxis([-0.2,1]);
end

fprintf('%10s %12.4e\n','qmin',qmin);
fprintf('%10s %12.4e\n','qmax',qmax);


if (Frame == 1 && mq == 3)
    cmax = max(abs([qmin,qmax]));
    caxis([-cmax,cmax]);
end    

hold on;
plot_stars();
hold off;



% cm = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0];
% colormap(cm);

colormap(parula);

colorbar;

NoQuery = 0;
prt = false;
if (prt)
    hidepatchborders;
    if (ishandle(circle))
        delete(circle);
    end
    mx = 8;
    maxlevel = 7;
    dpi = 2^7;    % fix at 128
    
    eff_res = mx*2^maxlevel;
    figsize = (eff_res/dpi)*[1,1];
    prefix = 'plot';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
end

shg

clear afterframe
