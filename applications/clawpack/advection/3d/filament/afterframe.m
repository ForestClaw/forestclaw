

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

yrbcolormap
caxis([0,1])
colorbar;

parms = read_vars();

% Compute exact solution at time T.

plot_soln = true;
if plot_soln
    N = 500;
    if (t > 0)
        [xout,yout] = filament_soln(N,t);
    else
        R_init = 0.25;
        xc_init = 0.5;
        yc_init = 1.0;
        th = linspace(0,2*pi,N+1);
        xout = R_init*cos(th) + xc_init;
        yout = R_init*sin(th) + yc_init;
    end
    hold on;
    for i = 1:length(zSliceCoords)  
        s = maxelev*zSliceCoords(i);
    
        zout = 0*xout + s;
        soln_handle = plot3(xout,yout,zout,'k','linewidth',3);
        % fprintf('Area of filament %24.16f\n',polyarea(xout,yout));
    end
    hold off;
end

showpatchborders;
setpatchborderprops('linewidth',1);
setpatchbordercolor('k');


view_3d = true;
if view_3d
    showcubes(5);
    showslices;
    view(3);
    daspect([1 1 1]);
    set(gca,'box','on');
    parms = read_vars();
    maxelev = parms.maxelev;
    axis([0 2 0 2 0 maxelev])
else
    daspect([1 1 1]);
    axis([0 2 0 2])
    view(2)
end

%axis off;
axis on;

shg;
