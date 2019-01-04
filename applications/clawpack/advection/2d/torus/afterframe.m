setviews;

alpha = 0.4;
s = 1e-2;
alim = [-1-alpha,1+alpha];
alim = alim + [-s s];
axis([alim alim]);

showpatchborders;
setpatchborderprops('linewidth',1);   % Default is 2
view(3);
daspect([1,1,1]);
axis off

isflat = true;
issquare = true;

yrbcolormap;
showpatchborders(1:10);
caxis([0,1])
qlo = 0;
qhi = 1;
under_label = sprintf('0 - %7.1e',qlo-qmin);
over_label = sprintf('1 + %7.1e',qmax-qhi);
fprintf('%6s %12s\n','qmin',under_label);
fprintf('%6s %12s\n\n','qmax',over_label);

hold on;
R = 0.5;  % revs per second
alpha = 0.4;
period = 0;
meqn = size(amrdata(1).data,1);
if (meqn == 1 && ~isflat)
    % assume we are doing the Gaussian problem with error
    [xp,yp,zp] = torus_soln(t,alpha,R,period);
    plot3(xp,yp,zp,'k','linewidth',2);
    hold off;
    hidepatchborders;
    view(3)
    camlight;
    axis off
end



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

showgridlines(1)
showpatchborders;
caxis([0,1e-10]);
view(2)

shg

clear afterframe;
clear mapc2m;
