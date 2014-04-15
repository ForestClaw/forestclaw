global cm_index;

if (mq == 1)
  cm_index = true;
else
  cm_index = false;
end

s = 1e-2;
axis([-s 2+s -s 2+s])
daspect([1 1 1]);
axis off;

if (mq == 1)
  hidegridlines;
  showpatchborders(1:6);

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
  yrbcolormap;
  showpatchborders(1:6);
  % hidepatchborders;
  delete(get(gca,'title'));
  caxis([0 1]);  
end;


view(2);

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'swirl0000','png');
  print('-dpng',filename);
end;

clear afterframe;
