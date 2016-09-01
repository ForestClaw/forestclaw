% bowl parameters
zmin = 80;
ep = 0.01;

if (PlotType == 1)
    s = 1e-2;
    axis([-100 100 -100 100])
    daspect([1 1 1]);
    axis on;
    
    showpatchborders(1:6);
    hidepatchborders(7)
    cv = linspace(qmin,qmax,20);
    % drawcontourlines(cv);
    set(gca,'zlim',[-8,0]);   % Need so that all patchborders show up
    
    rybcolormap;
    tol = 5;
    caxis([-tol,tol])
    % caxis([qmin,qmax])
    daspect([20,20,1]);
    
    hold on;
    zmin = 80;
    ep = 0.01;
    w = sqrt(zmin/ep);
    th = linspace(0,2*pi,500);
    plot(w*cos(th),w*sin(th),'k','linewidth',2);
    hold off;
    view(2);
    colorbar
else
    hold on;
    r = linspace(0,100,500);
    r2 = r.^2;
    b = ep*r.^2 - zmin;
    b(b >= 0) = 0;
    plot(r,b,'k','linewidth',2);  
    xlim([0,100]);
    ylim([-zmin,20]);    
    daspect([1 1 1])    
    [h_amr, labels_amr] = getlegendinfo;
    h = legend(h_amr,labels_amr,'location','southeast');
    hold off;
end

NoQuery = 0;
prt = false;
if (prt)
    filename = framename(Frame,'bowl0000','png');
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
