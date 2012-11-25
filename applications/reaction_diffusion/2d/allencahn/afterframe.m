s = 1e-2;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis off;
yrbcolormap;
% colormap(jet);

cv = 0.1:0.1:0.9;
% drawcontourlines(cv);
%showgridlines([1:4]);
hidepatchborders(6:7);
showpatchborders();
setpatchborderprops(1:8g,'linewidth',2);
caxis([-1 1]);

% view([-40 20])
% set(gca,'zdir','reverse');
% daspect([1 1 2]);

fprintf('[%e, %e],\n',qmin,qmax);

view(2);

NoQuery = 0;
MaxFrames = 120;
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
