

s = 1e-2;
axis([-1-s 1+s -1-s 1+s])
daspect([1 1 1]);
axis off;

if (mq == 1)
  yrbcolormap;
  caxis([0 1]);
  cm_index = false;  % global variable set in setprob.m
else
  cm_index = true;
  npmax = 10;
  cm = rand(npmax,3);
  colormap(cm);
  colorbar;
  caxis([1 npmax+1]);
end

setpatchborderprops(1:7,'linewidth',2);
showpatchborders;

view(2);

NoQuery = 0;
prt = false;
if (prt)
  filename = 'swirl000.png';
  str = num2str(Frame);
  len = length(str);
  filename(8-len+1:8) = str;
  pstr = ['print -dpng ',filename];
  disp(pstr);
  eval(pstr);
end;

clear afterframe;
clear setcolors;
