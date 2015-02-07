if (PlotType == 1)
  daspect([1 1 1]);
  rybcolormap;
  showgridlines(1:2);
  showpatchborders;

  caxis([-2 2]);
  colorbar;
elseif (PlotType == 4)
  hold on;
  dir = './1drad/_output/';
  dim = 1;
  [amrdata1d,t1d] = readamrdata(dim,Frame,dir);
  if (abs(t1d - t) > 1e-5)
    error('afterframe : 1d reference solution is not time synchronized');
  end;
  [q1d,x1d,p] = plotframe1ez(amrdata1d,mq,'b-');
  set(p,'linewidth',2);

  [ph,level_str] = getlegendinfo(0);  % base_level = 0
  lh = legend([ph,p],{level_str{:},'Exact'});
  set(lh,'fontsize',16);
  ylim([-2 3]);

  hold off;
end

view(2);

prt = false;
if (prt)
  fname = framename(Frame,'radial0000','png');
  print('-dpng',fname);
end

shg;

clear afterframe;
