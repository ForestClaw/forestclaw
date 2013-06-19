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


axis image;
showpatchborders;
setpatchborderprops(1:8,'linewidth',2);
showgridlines(1:3);
colormap(winter);

setviews;
% u = rrot(:,1);
% v = rrot(:,2);
% w = rrot(:,3);
% setcam(20*v);

prt = false;
MaxFrames  = 30;
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
