example = read_vars();

s = 0.0;
axis([-s, 1+s, -s, 1+s])
axis square;

fprintf('%-10s %16.8e\n','qmin',qmin);
fprintf('%-10s %16.8e\n','qmax',qmax);

if (ShowUnderOverShoots == 1)
    qlo = 0;
    qhi = 1;
    under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
    over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
    colorbar_underover(under_label,over_label);
end
if (mq == 3)
    % Plot the error
    c = max([qmin,qmax]);
    c = 3e-3;
    ca = [-c,c];    
else    
    % Plot the solution
    ca = [-0.3,1];
end

colormap(parula);
if (mq == 3)
    %colorbar
end
caxis(ca);

showpatchborders;
setpatchborderprops('linewidth',1);
hidegridlines;

if (example == 0)
    hold on;
    refine_threshold = 0.05;
    u = 1;   % (u,v) = velocity
    v = 0.5;
    maxlevel = 5;  % Doesn't have to match true maxlevel
    plot_refine_contours(mx,maxlevel,t,u,v,refine_threshold);
    hold off;
end

view(2)

% This is used for creating vectorized PDFs
prt_tikz = true;
if (prt_tikz)
    hidepatchborders;
    figsize = [16,16];  % Should match tikz figsize.
    maxlevel = 6;
    dpi = mx*2^maxlevel/figsize(1);
    prefix = 'plot';    
    caxis(ca);
    plot_tikz_fig(Frame,figsize,prefix,dpi);    
end



NoQuery = 0;
prt = true;
if (prt)
    filename = framename(Frame,'square0000','png');
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
