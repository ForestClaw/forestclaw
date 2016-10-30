axis image;
daspect([1 1 1]);
axis off;

yrbcolormap;
showpatchborders;
setpatchborderprops(1:8,'linewidth',1);
hidepatchborders(7);
% hidepatchborders;
view([65.5,12]);

% showgridlines(3);

setviews;
NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 41;
  filename = framename(Frame,'ccb0000','png');
  fprintf('Printing file %s\n',filename);
  print('-dpng','-r800',filename);
end

clear afterframe;
clear mapc2m;
clear mapc2m_pillow;
