% bowl parameters
zmin = 80;
ep = 0.01;

if (PlotType == 1)
    s = 1e-2;
    axis([-100 100 -100 100])
    daspect([1 1 1]);
    axis on;
    
    showpatchborders(1:9);
    cv = linspace(qmin,qmax,20);
    % drawcontourlines(cv);
    set(gca,'zlim',[-20,0]);   % Need so that all patchborders show up
    
    rybcolormap;
    tol = 1;
    caxis([-tol,tol])
    % caxis([qmin,qmax])
    daspect([20,20,1]);
    
    hold on;
    zmin = 80;
    ep = 0.01;
    w = sqrt(zmin/ep);
    th = linspace(0,2*pi,500);
    plot(w*cos(th),w*sin(th),'k','linewidth',2);
    % axis([50,75,50,75])
    % axis([0 100 0 100]);
    
    % Plot regions
    clear hh
    hh = plot_region([0,0],[100,100]);
    set(hh,'linewidth',3,'color','b');
    
    clear hh;
    if (t >= 3)
        hh(1) = plot_region([52,52],[72,72]);
        hh(2) = plot_region([75,-10],[95,10]);
        set(hh,'linewidth',3,'color','r');
    end
    
    clear hh;
    if (t >= 3.4)
        hh(1) = plot_region([57,57],[68,68]);
        hh(2) = plot_region([83,-4],[92,4]);
        set(hh,'linewidth',3,'color','w')
    end
    
    hold off;
    colorbar
elseif (PlotType == 4)
    % Plot 1d solution
    hold on;
    dir = './1drad/';
    dim = 1;
    [amrdata1d,t1d] = readamrdata(dim,Frame,dir);
    userfile = 'bathy1d';  % UserVariableFile
    [q1d,x1d] = plotframe1ez(amrdata1d,mq,'k-',userfile);

    % add bathymetry
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
