setviews;

global map isflat;

alpha = 0.4;
s = 1e-2;    
alim = [-1-alpha,1+alpha];
axis([alim alim]);
daspect([1 1 1]);
view(vtop)

yrbcolormap;

showpatchborders(1:10);
setpatchborderprops('linewidth',1)
caxis([0,1])


hidepatchborders(9)
showpatchborders;


%%
NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 8;
  axis([0 1 0 1]);
  filename = sprintf('annulus_%04d.png',Frame)
  print('-dpng',filename);
end

clear afterframe;
clear mapc2m;
