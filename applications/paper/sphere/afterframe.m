setviews;

global flat

example = read_vars();

if (flat)
    axis([0, 2*pi,  -pi/2, pi/2]);
else
    s = 0.0;
    axis([-1-s 1+s -1-s 1+s])
end    
daspect([1 1 1]);
axis on;

fprintf('%-10s %16.8e\n','qmin',qmin);
fprintf('%-10s %16.8e\n','qmax',qmax);

if (ShowUnderOverShoots == 1)
    qlo = 0;
    qhi = 1;
    under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
    over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
    colorbar_underover(under_label,over_label);
end
if (mq == 3)
    % Plot the error
    c = max([qmin,qmax]);
    ca = [-c,c];    
    ca = [-1,1]*1e-3;
else
    if (example == 3)
        if Frame > 0 && Frame < 10
            ca = [0,3.1];
        else
            ca = [-0.2,10];
        end
    else
        ca = [0,1];
    end
end

colormap(parula);
if (mq == 3)
    %colorbar
end
caxis(ca);

showpatchborders;
setpatchborderprops('linewidth',1);
hidegridlines;

if (~flat)
    view([ 59.297571594931483, 8.368467153284623]);
    view([1.527452862277914e+02, 1.780446485025879]);
    view(3)
    % view(vtop);
    % view(vtop);
    hold on;
    th = linspace(-pi/2, pi/2,500);
    plot3(cos(th),0*th,sin(th),'k','linewidth',4);
    th = linspace(pi/32, 2*pi-pi/32,500);
    plot3(cos(th),sin(th),0*th,'k','linewidth',4);
    hold off;
end

cv = linspace(0.1,3.1,20);
cv([1 end]) = [];
% drawcontourlines(cv);

NoQuery = 0;
prt = false;
if (prt)
    filename = framename(Frame,'swirl0000','png');
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
