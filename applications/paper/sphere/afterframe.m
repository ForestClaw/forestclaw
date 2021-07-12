setviews;

global flat

example = read_vars();

% Data info
fprintf('%-10s %16.8e\n','qmin',qmin);
fprintf('%-10s %16.8e\n','qmax',qmax);

if (ShowUnderOverShoots == 1)
    qlo = 0;
    qhi = 1;
    under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
    over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
    colorbar_underover(under_label,over_labcel);
end
if (mq == 3)
    % Plot the error
    c = max([qmin,qmax]);
    ca = [-1,1]*1e-4;
else
    if (example == 3)
        if Frame > 0 && Frame < 10
            ca = [0,3.1];
        else
            ca = [-0.2,10];
        end
    else
        ca = [-0.2,1];
    end
end

% Color axis
colormap(parula);
if (mq == 3)
    colorbar
end
caxis(ca);

% AMR Patches
showpatchborders;
setpatchborderprops('linewidth',1);
hidegridlines;

% Grid axes
if (flat)
    axis([0, 2*pi,  -pi/2, pi/2]);
else
    s = 0.0;
    axis([-1-s 1+s -1-s 1+s])
end    
daspect([1 1 1]);
axis on;
set(gca,'fontsize',16);


if (~flat)
    view([ 59.297571594931483, 8.368467153284623]);
    view([1.527452862277914e+02, 1.780446485025879]);
    view(3)
    % view(vtop);
    % view(vtop);
    hold on;
    th = linspace(-pi/2, pi/2,500);
    % plot3(cos(th),0*th,sin(th),'k','linewidth',4);
    th = linspace(pi/32, 2*pi-pi/32,500);
    % plot3(cos(th),sin(th),0*th,'k','linewidth',4);
    hold off;
end

cv = linspace(0.1,3.1,20);
cv([1 end]) = [];
% drawcontourlines(cv);

%{
delete(get(gca,'title'));
o = findobj('Type','colorbar');
delete(o);
set(gca,'clipping','off')
view([-37.5000, 30]);
zoom(1.5);
axis off;
%}

view(2);

NoQuery = 0;
MaxFrames = 10;
prt = true;
if (prt)
    filename = sprintf('maxlevel3/errs_maxlevel3_sync_Frame%02d',Frame);
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
