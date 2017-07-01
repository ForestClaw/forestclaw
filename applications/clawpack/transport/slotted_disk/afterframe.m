axis image;
daspect([1 1 1]);
axis off;

yrbcolormap;

% showgridlines(1:4);
showpatchborders;
setpatchborderprops('linewidth',1);
hidepatchborders(7);
view([65.5,12]);

setviews;
% set(gcf,'visible','off');

NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 10;  
  filename = framename(Frame,'sdisk000','png');
  fprintf('Printing file %s\n',filename);
  print('-dpng','-r600',filename);
end

shg

clear afterframe;
clear mapc2m;
