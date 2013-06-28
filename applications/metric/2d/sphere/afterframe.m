s = 1e-2;
% axis([-1-s 1+s -1-s 1+s -1-s 1+s])
axis image;
daspect([1 1 1]);
axis off;

colormap(jet);
caxis([qmin qmax]);

showgridlines(1:3);
setpatchborderprops(1:6,'linewidth',2);
showpatchborders;
setviews;

colorbar;

NoQuery = 0;
prt = false;
if (prt)
  filename = 'sphere000.png';
  str = num2str(Frame);
  len = length(str);
  filename(8-len+1:8) = str;
  pstr = ['print -dpng ',filename];
  disp(pstr);
  eval(pstr);
end;

clear afterframe;
clear mapc2m;
clear mapc2m_pillow;
