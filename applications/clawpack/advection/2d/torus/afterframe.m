setviews;

alpha = 0.4;
s = 1e-2;
alim = [-1-alpha,1+alpha];
alim = alim + [-s s];
axis([alim alim]);
daspect([1 1 1]);

yrbcolormap;

showpatchborders;
setpatchborderprops('linewidth',1)
view(3);
axis off

shg

clear afterframe;
clear mapc2m;
