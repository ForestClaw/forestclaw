if (PlotType == 1)
  rybcolormap;
  caxis([-2 2]);
  colorbar;  

  showpatchborders;
  setpatchborderprops('linewidth',1);

  axis off;
  daspect([1 1 1]);
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

  plot([1 1],ylim,'k--');
  
  hold off;
end

view(2);
set(gca,'fontsize',16)

prt = false;
if (prt)
  fname = framename(Frame,'radial0000','png');
  print('-dpng',fname);
end

shg;

clear afterframe;
