s = 0;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis on;

ex = 2;
lstr = 'no limiter';
mstr = 'QS';


if (ShowUnderOverShoots == 1)
    qlo = 0;
    qhi = 1;
    under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
    over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
    fprintf('%-10s %12.4e\n','qmin',qmin);
    fprintf('%-10s %12.4e\n','qmax',qmax);
    colorbar_underover(under_label,over_label);
else
    if (ex == 1)
        ca = [0 9];
    else
        ca = [0 16];
    end    
    colormap(parula);
    colorbar;
    caxis(ca);
    % title(sprintf('Example %d',ex),'fontsize',18);
    title(sprintf('%s (%s)',mstr,limstr),'fontsize',18);

    
    figure(2);
    clf;
    [xm,ym] = meshgrid(xcenter,ycenter);
    h = surf(xm,ym,q,'edgecolor','none');
    caxis(ca);
    set(gca,'xlim',[0 1],'ylim',[0 1],'zlim',ca,'box','on');    
    colorbar;
    camlight;
    xlabel('x','fontsize',16);
    ylabel('y','fontsize',16);
    % title(sprintf('Example %d',ex),'fontsize',18);
    title(sprintf('%s (%s)',mstr,limstr),'fontsize',18);
    
    figure(1);    
end
fprintf('%-10s %12.4e\n','qmin',qmin);
fprintf('%-10s %12.4e\n','qmax',qmax);


showpatchborders;
setpatchborderprops('linewidth',1);
view(2);

set(gca,'fontsize',16);
set(gca,'box','on');

NoQuery = 0;
prt = false;
if (prt)
    filename = framename(Frame,'swirl0000','png');
    print('-dpng',filename);
end

shg

clear afterframe;
