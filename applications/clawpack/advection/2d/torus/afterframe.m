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

%{
view(2)
axis([0,1,0,1]);
xtick = linspace(0,1,9);
set(gca,'xtick',xtick);
set(gca,'ytick',xtick);
set(gca,'xticklabels',2*(xtick-0.5));
set(gca,'yticklabels',(xtick-0.5));
daspect([2 5 1]);
axis on
%}

% showgridlines(1:3)

shg

clear afterframe;
clear mapc2m;
