s = 0;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis on;

if (ShowUnderOverShoots == 1)
    qlo = 0;
    qhi = 1;
    under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
    over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
    fprintf('%-10s %12.4e\n','qmin',qmin);
    fprintf('%-10s %12.4e\n','qmax',qmax);
    colorbar_underover(under_label,over_label);
else
    colormap(parula);
    colorbar;
    caxis([0 5]);
    
    figure(2);
    clf;
    h = surf(q,'edgecolor','none');
    camlight;
    caxis([0 5]);
    colorbar;
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
