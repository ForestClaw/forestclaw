axis image;
daspect([1 1 1]);
axis off;

yrbcolormap;
qlo = 0.1;
qhi = 1.1;
% showgridlines(1:4);
showpatchborders(1:8);
setpatchborderprops(1:8,'linewidth',1);
view([-51.5,0]);

setviews;
NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 21;
  filename = framename(Frame,'ccb000','png');
  fprintf('Printing file %s\n',filename);
  print('-dpng','-r800',filename);
end;

clear afterframe;
clear mapc2m;
clear mapc2m_pillow;
