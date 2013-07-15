s = 1e-2;
axis([-1-s 1+s -1-s 1+s])
daspect([1 1 1]);
axis off;

if (mq == 1)
    yrbcolormap;
    caxis([0 1]);
  else
    cm = colorcube;
    colormap(cm(1:end-10,:));
    colorbar;
    caxis([qmin qmax]+1);
end;

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
