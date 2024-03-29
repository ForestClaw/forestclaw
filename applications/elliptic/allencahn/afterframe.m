
% Solution or error data
fprintf('%10s %12.4e\n','qmin',qmin);
fprintf('%10s %12.4e\n','qmax',qmax);

% Mesh features
showpatchborders;
setpatchborderprops('linewidth',1);

ca = [-1.5,1];
caxis(ca);

% Contour lines
cv = linspace(0,1,11);
% drawcontourlines(cv);


% Color map and axis
colormap(parula);
colorbar;

% Axis labels
ax = 0;
bx = 1;
ay = 0;
by = 1;
s = 1e-2;
axis([ax-s bx+s ay-s by+s])
daspect([1 1 1]);

xlabel('x','fontsize',16);
ylabel('y','fontsize',16);


% Surface plot
if (length(amrdata) == 1)
    % We are only on a single grid
    figure(2);
    clf;    
    h = surf(xcenter,ycenter,q');
    set(h,'edgecolor','none');
    view(3);
    axis square;    
    % view(vfront);
    set(gca,'zlim',[-1.1,1.1]);
%     set(gcf,'color','k');
%     set(gcf,'clipping','off');
%     axis off;
%     camlight;
    figure(1);
end


NoQuery = 0;
prt = false;
MaxFrames = 10000;
if (prt)
    % hidepatchborders;
    fname = sprintf('allencahn%04d.png',Frame);
    print(fname,'-dpng');
    fprintf('Printing %s\n',fname);
%     if (ishandle(circle))
%         delete(circle);
%     end
%     mx = 16;
%     maxlevel = 6;
%     dpi = 2^7;    % fix at 128
%     
%     eff_res = mx*2^maxlevel;
%     figsize = (eff_res/dpi)*[1,1];
%     prefix = 'plot';
%     plot_tikz_fig(Frame,figsize,prefix,dpi)
end




shg

clear afterframe
