setviews;

example = read_vars();

s = 0.0;
axis([-1-s 1+s -1-s 1+s])
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
else    
    % Plot the solution
    ca = [-0.3,1];
end

colormap(parula);
if (mq == 3)
    %colorbar
end
caxis(ca);

showpatchborders;
setpatchborderprops('linewidth',1);
hidegridlines;

view([ 59.297571594931483, 8.368467153284623]);
view([1.527452862277914e+02, 1.780446485025879]);
view(3)
view(vleft);
% view(vtop);

NoQuery = 0;
prt = false;
if (prt)
    filename = framename(Frame,'swirl0000','png');
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
