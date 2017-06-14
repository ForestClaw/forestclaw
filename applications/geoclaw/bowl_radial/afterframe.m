% bowl parameters
zmin = 80;
ep = 0.01;

if (PlotType == 1)
    showpatchborders;
    setpatchborderprops('linewidth',2);
    cv = linspace(qmin,qmax,20);
    % drawcontourlines(cv);
    set(gca,'zlim',[-20,0]);   % Need so that all patchborders show up
    
    colormap(parula);
    tol = 1;
    caxis([-tol,tol])
    daspect([20,20,1]);
    colorbar
    
    hold on;
    zmin = 80;
    ep = 0.01;
    w = sqrt(zmin/ep);
    th = linspace(0,2*pi,500);
    plot(w*cos(th),w*sin(th),'k','linewidth',2);
    
    add_gauges();
    add_regions(t);
    
    hold off;
    
    s = 1e-2;
    axis([-100 100 -100 100])
    daspect([1 1 1]);
    set(gca,'fontsize',16);
    axis on;
    
elseif (PlotType == 4)
    % Plot 1d solution
    hold on;
    dir = './1drad_pyclaw/_output';
    dim = 1;
    [auxdata1d,] = readamrdata(dim,0,dir,'aux');
    [amrdata1d,t1d] = readamrdata(dim,Frame,dir);
    userfile = 'bathy1d';  % UserVariableFile
    [q1d,x1d] = plotframe1ez(amrdata1d,mq,'k-','','',auxdata1d);

    % add bathymetry
    r = linspace(0,100,500);
    r2 = r.^2;
    b = ep*r.^2 - zmin;
    b(b >= 0) = 0;
    plot(r,b,'k','linewidth',2);  
    xlim([0,100]);
    ylim([-80,40]);
    daspect([10 10 1])
    [h_amr, labels_amr] = getlegendinfo(0);   % Input base level
    h = legend(h_amr,labels_amr,'location','southwest');
    set(h,'fontsize',16);
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
