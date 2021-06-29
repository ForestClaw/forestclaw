setviews;

yrbcolormap;
caxis([0 1]);

showpatchborders;
setpatchborderprops('linewidth',1);

daspect([1 1 1]);
axis off;
view([117.50,4]);
% caxis([0.75,1.5]);

NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 41;
  filename = framename(Frame,'gaussian0000','png');
  fprintf('Printing file %s\n',filename);
  print('-dpng','-r800',filename);
end

clear afterframe;
clear mapc2m;
