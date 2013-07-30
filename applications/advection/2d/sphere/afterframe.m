global cm_index;
s = 1e-2;
axis image;
daspect([1 1 1]);
axis off;


if (mq == 1)
  cm_index = true;
  hidegridlines;

  npmax = 10;
  if (Frame == 0)
    cm = multicolormap(npmax);
  end;
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

  setpatchborderprops(1:10,'linewidth',2);
  showpatchborders;

else
  cm_index = false;

  yrbcolormap;
  caxis([0 1]);

  showgridlines(1:2);
  setpatchborderprops(1:6,'linewidth',2);
  hidepatchborders;
  showpatchborders(1:5);
end
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

clear afterframe;
clear mapc2m;
clear mapc2m_pillow;
