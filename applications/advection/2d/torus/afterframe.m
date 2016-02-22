setviews;

global map isflat;

alpha = 0.4;
s = 1e-2;    
if strcmp(map,'nomap')
    axis([-1 1 -1 1]);
else
    if strcmp(map,'cart')
        alim = [-1 1];
    elseif strcmp(map,'latlong')
        alim = [-1 1];
    else
        alim = [-1-alpha,1+alpha];
    end
    alim = alim + [-s s];
    axis([alim alim]);
    daspect([1 1 1]);
end

yrbcolormap;
showpatchborders(1:10);
caxis([0,1])
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
daspect([1,1,1]);

if (isflat)
    view(2);
else
    hidepatchborders;
    view(3);
    % camlight;
    showpatchborders;
    hold on;
    R = 0.5;  % revs per second
    alpha = 0.4;
    period = 16;
    [xp,yp,zp] = torus_soln(t,alpha,R,period);    % camlight;
    plot3(xp,yp,zp,'k','linewidth',2);
    hold off;
end

setpatchborderprops('linewidth',1)
hidepatchborders(9)

NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 128;
  filename = framename(Frame,'torus0000','png');
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
clear mapc2m_torus;
clear torus_soln;
clear parallelpartitions;
