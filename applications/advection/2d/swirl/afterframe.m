global cm_index;

global cm_index;

s = 1e-2;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis off;

if (mq == 1)
  hidegridlines;

  npmax = 10;
  if (Frame == 0)
    cm = multicolormap(npmax);
  end;
  colormap(cm);
  caxis([1 npmax+1]);
  showpatchborders;

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
  showpatchborders(1:10);
  caxis([0,1])
  qlo = 0;
  qhi = 1;
  under_label = sprintf('0 - %7.1e',qlo-qmin);
  over_label = sprintf('1 + %7.1e',qmax-qhi);
  fprintf('%6s %12s\n','qmin',under_label);
  fprintf('%6s %12s\n\n','qmax',over_label);

  if (ShowUnderOverShoots)
    qlo = 0;
    qhi = 1;
    colorbar_underover(under_label,over_label);
  end
end


view(2);

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'swirl0000','png');
  print('-dpng',filename);
end

shg
clear afterframe;
