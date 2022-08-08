yrbcolormap;
setviews;

showpatchborders(1:10);
setpatchborderprops('linewidth',1)
caxis([0,1])

hideslices();
showslices('z',3);


parms = read_vars();
alpha = parms.maxelev + 0.005;
alpha = 0.5;
s = 1e-2;    
alim = [-1-alpha,1+alpha];
axis([alim alim alim]);

daspect([1 1 1]);
view(vtop);
%view(3);

shg
