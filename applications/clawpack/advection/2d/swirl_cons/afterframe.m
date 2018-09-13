s = 0.01;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis on;

ex = 1;
lstr = 'no limiter';
mstr = 'WD';

if (PlotType == 4)
    hold on;
    [grids, soln] = test_cons();
    mg =length(grids);
    solnq = soln{Frame+1}.q;
    for m = 1:mg
        g = grids{m};
        q = solnq{m};
        plot(g.xc(:),q(:),'bx','linewidth',2);
    end
    axis([0 1 0 3.1]);
    hold off;    
elseif (ShowUnderOverShoots == 1)
    qlo = 0;
    qhi = 1;
    under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
    over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
    fprintf('%-10s %16.8e\n','qmin',qmin);
    fprintf('%-10s %16.8e\n','qmax',qmax);
    colorbar_underover(under_label,over_label);
else
    switch ex
        case 1
            ca = [-1.4, 2];
        case 2
            ca = [-1, 2.5];    
        case 3
            ca = 1*[-1,1];
            % ca = 1e-15*[-1,1];
    end    
    %yrbcolormap;
    colormap(parula);
    colorbar;
    caxis(ca);
    % title(sprintf('Example %d',ex),'fontsize',18);
    title(sprintf('%s (%s)',mstr,lstr),'fontsize',18);

    cv = linspace(ca(1),ca(2),11);
    % drawcontourlines(cv);
    
    one_patch = length(amrdata) == 1;
    if (one_patch)
        figure(2);  
        clf;    
        [xm,ym]     = meshgrid(xcenter,ycenter);
        h = surf(xm,ym,q,'edgecolor','none');
        caxis(ca);
        set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',ca,'box','on');
        colorbar;
        camlight;
        xlabel('x','fontsize',16);
        ylabel('y','fontsize',16);
        % title(sprintf('Example %d',ex),'fontsize',18);
        title(sprintf('%s (%s)',mstr,lstr),'fontsize',18);
        figure(1);    
        
    end
    setpatchborderprops('linewidth',1);
    showpatchborders;
    daspect([1 1 1]);
    view(2);        
    
    fprintf('%-10s %16.8e\n','qmin',qmin);
    fprintf('%-10s %16.8e\n','qmax',qmax);    
end

% showgridlines
% hidegridlines(5);
% hidepatchborders(6)
% setpatchborderprops('linewidth',2);
hidepatchborders;
hidegridlines;

showgridlines;
showpatchborders;




set(gca,'fontsize',16);
% set(gca,'box','on');

NoQuery = 0;
prt = false;
if (prt)
    filename = framename(Frame,'swirl0000','png');
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
