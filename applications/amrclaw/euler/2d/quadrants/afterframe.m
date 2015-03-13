s = 1e-2;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis off;

colormap(jet);
cv = linspace(qmin,qmax,21);
cv([1 end]) = [];
% drawcontourlines(cv);
% setcontourlineprops('linewidth',2);
setpatchborderprops('linewidth',1);
showpatchborders;
colorbar;


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

shg

clear afterframe;
