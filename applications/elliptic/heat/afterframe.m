
% Solution or error data
fprintf('%10s %12.4e\n','qmin',qmin);
fprintf('%10s %12.4e\n','qmax',qmax);

% Mesh features
showpatchborders;
setpatchborderprops('linewidth',1);

example = read_vars();
switch example
    case 1
        ca = [0.9,2];
    case 4
        ca = [0,1];
end
ca = [qmin,qmax];
caxis(ca);
cv = linspace(ca(1),ca(2),11);
drawcontourlines(cv);

% Color map and axis
colormap(parula);
colorbar;

% Add some extra info
% hold on;
% plot_stars();
% hold off;

% Surface plot
if (length(amrdata) == 1)
    % We are only on a single grid
    figure(2);
    clf;    
    h = surf(xcenter,ycenter,q');
    set(h,'edgecolor','none');
    view(3);
    axis square;
    % set(gca,'zlim',ca);
%     set(gcf,'color','k');
%     set(gcf,'clipping','off');
%     axis off;
%     camlight;
    figure(1);
end

% Axis labels
ax = -1;
bx = 1;
ay = -1;
by = 1;
s = 1e-2;
axis([ax-s bx+s ay-s by+s])
daspect([1 1 1]);

xlabel('x','fontsize',16);
ylabel('y','fontsize',16);

% title(tstr,'fontsize',18);

        
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
