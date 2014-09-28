axis image;
daspect([1 1 1]);

colormap(jet);

showpatchborders;
setviews;

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

shg;

clear afterframe;
clear mapc2m;
clear mapc2m_pillow;
