setviews;

alpha = 0.4;
period = 0;

if (ShowUnderOverShoots)
    qlo = 0;
    qhi = 1;
    colorbar_underover(under_label,over_label);
end

yrbcolormap;
showpatchborders;
setpatchborderprops('linewidth',1);
daspect([1,1,1]);

caxis([0 1])

o = findobj('type','Light');
if ishandle(o)
    delete(o);
end
%camlight

view(3);
% view(vtop);

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
