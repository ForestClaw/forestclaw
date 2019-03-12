s = 0.01;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis on;

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
    ca = [-max([qmin,qmax]),max([qmin,qmax])];
else    
    % Plot the solution
    ca = [-1.4, 2];
end

colormap(jet);
% colorbar;
caxis(ca);

showpatchborders;
setpatchborderprops('linewidth',1);
hidegridlines;

hold on;
refine_threshold = 0.3;
u = 1;
v = 1;
maxlevel = 5;  % Doesn't have to match true maxlevel
plot_refine_contours(mx,maxlevel,t,u,v,refine_threshold);
hold off;

% This is used for creating vectorized PDFs
prt_tikz = false;
if (prt_tikz)
    figsize = [64,64];  % Should match tikz figsize.
    maxlevel = 8;
    dpi = mx*2^maxlevel/figsize(1);
    prefix = 'plot';    
    % caxis([-1,1]*5e-6)
    caxis(ca/10);
    plot_tikz_fig(Frame,figsize,prefix,dpi);    
end


NoQuery = 0;
prt = false;
if (prt)
    filename = framename(Frame,'swirl0000','png');
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
