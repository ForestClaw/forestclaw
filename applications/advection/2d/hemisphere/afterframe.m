s = 1e-2;
% axis([-1-s 1+s -1-s 1+s -1-s 1+s])
axis image;
daspect([1 1 1]);
axis off;
yrbcolormap;
% colormap(jet);

cv = 0.1:0.1:0.9;
% drawcontourlines(cv);
caxis([0 1]);

showgridlines(1:3)
setpatchborderprops(1:7,'linewidth',2);
% hidepatchborders(7);

% hidepatchborders;
% hidegridlines;

% view([-121.5 30]);
setviews;
view(vtop);
view(3);

s = sum(sum(isnan(q)));
fprintf('Number of nans : %d\n',s);

NoQuery = 0;
MaxFrames = 192;
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
clear mapc2m;
