draw_mesh = false;
draw_scatter = false;

if (draw_scatter)
  view(2);
  set(gca,'ylim',[-2 2]);
  set(gca,'xlim',[-pi/2 pi/2]);
  clear map1d;
  clear afterframe;
  return;
end

axis(rsphere*1.1*[-1 1 -1 1 -1 1]);
daspect([1 1 1]);
axis off;

if (~draw_mesh)
  c1 = -1;
  c2 = 3;
  s = 0.001;
  ca = [c1-s c2+s];
  colormap(winter);

  caxis(ca);

  cv = linspace(ca(1),ca(2),11);
  cv(1) = [];
  cv(2) = [];
  drawcontourlines(cv);
  showpatchborders;
  setpatchborderprops(1:6,'linewidth',2);
  % showgridlines(1:4);
else
  colormap(white);
  showpatchborders;
  showgridlines([1 2]);
end;

setviews;

prt = false;
if (prt)
  NoQuery = false;
  nstr = num2str(Frame);
  len = length(nstr);
  if (draw_mesh)
    fstr = 'swe_mesh0000.png';
    fstr(12-len+1:12) = nstr;
  else
    fstr = 'swe0000.png';
    fstr(7-len+1:7) = nstr;
  end;
  pstr = ['print -dpng ',fstr];
  disp(pstr);
  eval(pstr);
else
  NoQuery = false;
end;


clear afterframe;
clear mapc2m;
