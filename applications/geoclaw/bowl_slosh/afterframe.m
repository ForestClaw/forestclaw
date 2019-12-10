% bowl parameters
zmin = 80;
ep = 0.01;

if (PlotType == 1)
    hold on;
    
    % Patches
    showgridlines(1:6);
    showpatchborders(1:6);
    % hidepatchborders(5);
    setpatchborderprops('linewidth',2);
    set(gca,'zlim',[-20 1]);   % Need so that all patchborders show up
    
    colormap(parula);
    colorbar;
    tol = -0.8;
    c1 = -0.14;
    c2 = 0.15;
    caxis([c1,c2]);
    
    % Contour lines (some are hidden by patch boundaries)
    cv = linspace(c1,c2,21);
    drawcontourlines(cv);
    
    % Plot boundary of true solution
    a = 1;
    sigma = 0.5;
    h0 = 0.1;
    grav = 9.81;
    omega = sqrt(2*grav*h0) / a;

    xe = linspace(-2,2,200);
    ye = linspace(-2,2,200);
    [xem,yem] = meshgrid(xe,ye);
    
    y = 0;
    B = h0*(xem.^2 + yem.^2)/a^2 - h0;
    eta1 = sigma*h0/a.^2 * (2*xem*cos(omega*t) + 2.*yem*sin(omega*t) -sigma);
    
    % contour(xem,yem,B-eta1,[0 0],'k','linewidth',5);
   
    hold off;

    disp(qmin)
    disp(qmax)
    
    
    % Axes
    axis([-2 2 -2 2])
    daspect([1 1 1]);
    set(gca,'fontsize',16);
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

title(sprintf('slosh (fclaw) at time %g',t),'fontsize',18);

NoQuery = 0;
prt = false;
MaxFrames = 16;
if (prt)
    filename = sprintf('bowl%04d.png',Frame);
    fprintf('Print file %s\n',filename);
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
