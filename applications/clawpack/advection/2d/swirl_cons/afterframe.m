s = 0.05;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis on;

ex = 3;
lstr = 'no limiter';
mstr = 'QS';

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
    fprintf('%-10s %12.4e\n','qmin',qmin);
    fprintf('%-10s %12.4e\n','qmax',qmax);
    colorbar_underover(under_label,over_label);
else
    switch ex
        case 0
            ca = [0 9];
        case 1
            ca = [0 40];
        case 2
            ca = [0 16];    
        case 3
            ca = [0,1.5];
%             ca = 1e-15*[-1,1];
    end    
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
    
    fprintf('%-10s %12.4e\n','qmin',qmin);
    fprintf('%-10s %12.4e\n','qmax',qmax);    
end



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
