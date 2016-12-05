setviews;

alpha = 0.4;
s = 1e-2;
alim = [-1-alpha,1+alpha];
alim = alim + [-s s];
axis([alim alim]);
daspect([1 1 1]);

if PlotParallelPartitions ~= 1
    yrbcolormap;
end

if (mq == 2)
    % plotting error
    qc = max(abs([qmin,qmax]));
    caxis([-qc,qc]);
    fprintf('%6s %12s\n\n','qmax',qmax);
    daspect([1,1,1]);
end

view(3);
% camlight;
showpatchborders;

hold on;
R = 0.5;  % revs per second
alpha = 0.4;
period = 16;
meqn = size(amrdata(1).data,1);
if (meqn == 2)
    % assume we are doing the Gaussian problem with error
    [xp,yp,zp] = torus_soln(t,alpha,R,period);
    plot3(xp,yp,zp,'k','linewidth',2);
    hold off;
    hidepatchborders;
    view(3)
    camlight;
    axis off
end

showpatchborders;
setpatchborderprops('linewidth',1)

axis off


shg

clear afterframe;
clear mapc2m;
clear mapc2m_torus;
clear torus_soln;
clear parallelpartitions;
