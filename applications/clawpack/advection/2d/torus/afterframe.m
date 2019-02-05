setviews;

global map_isflat;

isflat = map_isflat;
issquare = true;

plot_contour = false;
plot_qtrue = false;

R = 0.5;  % revs per second
alpha = 0.4;
period = 0;
meqn = size(amrdata(1).data,1);

if plot_qtrue
    evec = zeros(3,1);
    for i = 1:length(amrdata)
        ad = amrdata(i);
        xlow = ad.xlow;
        ylow = ad.ylow;
        xe = linspace(xlow,xlow+ad.mx*ad.dx,ad.mx+1);
        ye = linspace(ylow,ylow+ad.my*ad.dy,ad.my+1);
        xc = xe(1:end-1) + ad.dx/2;
        yc = ye(1:end-1) + ad.dy/2;
        q = reshape(ad.data(1,:),ad.mx,ad.my)';
        set_blocknumber(ad.blockno);
        [xe1,ye1,~] = mapc2m_brick(xe,ye,0);
        [xc1,yc1,~] = mapc2m_brick(xc,yc,0);
        qtrue = torus_soln_grid(xc1,yc1,xe1,ye1,t,alpha,R,period);  
        err = abs(qtrue(:)-q(:));
        evec(1) = evec(1) + sum(err)*ad.dx;
        evec(2) = evec(2) + sum(err.^2)*ad.dx;
        evec(3) = max([evec(3); err]);
        hold on;
    end
    area_exact = 4*alpha*pi^2;
    evec(1) = evec(1)/area_exact;
    evec(2) = sqrt(evec(2)/area_exact);
    fprintf('\n');
    fprintf('%5d %12.4e %12.4e %12.4e\n',mx,evec(1),evec(2),evec(3));
    fprintf('\n');
end

if plot_contour
    yrbcolormap;
    hold on;
    % assume we are doing the Gaussian problem with error
    [xp,yp,zp] = torus_soln_contour(t,alpha,R,period);
    plot3(xp,yp,zp,'k','linewidth',2);
    hold off;
    hidepatchborders;
    view(3)
    % camlight;
    axis off
end


qlo = 0;
qhi = 1;

under_label = sprintf('0 - %7.1e',qlo-qmin);
over_label = sprintf('1 + %7.1e',qmax-qhi);

fprintf('%6s %12s\n','qmin',under_label);
fprintf('%6s %12s\n\n','qmax',over_label);


if (ShowUnderOverShoots)
    qlo = 0;
    qhi = 1;
    colorbar_underover(under_label,over_label);
end

if isflat
    view(2)
    if (issquare)
        axis([0,1,0,1]);
        daspect([1,1,1]);
    else
        xtick = linspace(0,1,9);
        set(gca,'xtick',xtick);
        set(gca,'ytick',xtick);
        set(gca,'xticklabels',2*(xtick-0.5));
        set(gca,'yticklabels',(xtick-0.5));
        daspect([2 5 1]);
    end
    axis on
end

yrbcolormap;
% colormap(parula);
showgridlines(1:5);
showpatchborders;
% hidepatchborders;
setpatchborderprops('linewidth',1);
daspect([1,1,1]);

caxis([0 1])
% colorbar
axis off;

o = findobj('type','Light');
if ishandle(o)
    delete(o);
end
%camlight

% view(vfront);
view(3);

prt = false;
NoQuery = false;
if (prt)
    delete(get(gca,'title'));
    fname = sprintf('torus%02d.png',Frame);
    print('-dpng','-r640',fname);
end
shg

clear afterframe;
clear mapc2m;
