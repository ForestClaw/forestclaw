s = 1e-2;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis off;

if (mq == 1)
  % Plot partitions
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
else
  colormap(jet);
  cv = linspace(qmin,qmax,21);
  cv([1 end]) = [];
  drawcontourlines(cv);
  setcontourlineprops('linewidth',2);
  setpatchborderprops(1:7,'linewidth',2);
  showpatchborders;
  colorbar;
end;


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
