ax = 0;
bx = 1;
ay = 0;
by = 1;
s = 1e-2;
axis([ax-s bx+s ay-s by+s])
daspect([1 1 1]);
axis off;

showpatchborders;
setpatchborderprops('linewidth',1);

view(2);

rhs_choice = 0;
circle = [];
if rhs_choice == 3
    th = linspace(0,2*pi,500);
    r = 0.25;
    hold on;
    circle = plot(r*cos(th)+0.5,r*sin(th)+0.5,'k','linewidth',2);
    hold off;
end

% cm = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0];
% colormap(cm);

colormap(parula);

NoQuery = 0;
prt = true;
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
