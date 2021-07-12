yrbcolormap;
setviews;

showpatchborders(1:10);
setpatchborderprops('linewidth',1)
caxis([0,1])

alpha = 0.4;
s = 1e-2;    
alim = [-1-alpha,1+alpha];
axis([alim alim]);

daspect([1 1 1]);
view(3);

shg
