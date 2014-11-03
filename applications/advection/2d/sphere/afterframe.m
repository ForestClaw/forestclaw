setviews;
s = 1e-2;
daspect([1 1 1]);
axis off;


if (mq == 0)
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
  yrbcolormap;
  caxis([0 1]);
  showgridlines(1:2);
  showpatchborders(1:7);
  if (ShowUnderOverShoots == 1)
    under_label = sprintf('%6.2e',qmin);
    over_label = sprintf('1 + %6.2e',qmax-1);
    colorbar_underover(under_label,over_label);
  end
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
zoom(1.1);
zoom(1.1);
zoom(1.1);
zoom(1.1);
zoom(1.1);
zoom(1.1);
hidepatchborders;
hidegridlines;
showgridlines(1:6);
shg;


clear afterframe;
clear mapc2m;
clear mapc2m_cubed_sphere;
