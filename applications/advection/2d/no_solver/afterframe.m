

s = 1e-2;
axis([-1-s 1+s -1-s 1+s])
daspect([1 1 1]);
axis off;

if (mq == 1)
  yrbcolormap;
  caxis([0 1]);
  cm_index = false;  % global variable set in setprob.m
else
  cm_index = true;  % Used by setcolors.m
  npmax = 10;
  if (Frame == 0)
    cm = multicolormap(npmax);
  end
  colormap(cm);
  caxis([1 npmax+1]);

  % Fix colorbar
  colorbar;
  o = findobj('Tag','Colorbar');
  set(o,'ytick',(1:(npmax+1)) + 0.5);
  set(o,'yticklabel',(1:npmax)-1);
  set(o,'ylim',[qmin+1 qmax+2]);
  set(o,'ticklength',[0 0])
  set(o,'fontsize',16,'fontweight','bold')
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
