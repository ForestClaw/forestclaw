global boxes;

s = 1e-2;
% axis([-1-s 1+s -1-s 1+s -1-s 1+s])
axis image;
daspect([1 1 1]);
axis off;

colormap([1 1 0]);

showgridlines(1:3);
setpatchborderprops(1:8,'linewidth',1);
showpatchborders;
setviews;
view(2);


clear circle_h;
clear box_h;
for i = 1:length(boxes),
  [circle_h(i),diag_h,box_h(i)] = draw_circles(boxes(i));
  set(circle_h(i),'linewidth',3,'color','k');
  set(box_h(i),'linewidth',3,'visible','on');
  set(diag_h,'visible','off');
end;

% colorbar;

NoQuery = 0;
prt = false;
if (prt)
  filename = 'sphere000.png';
  str = num2str(Frame);
  len = length(str);
  filename(8-len+1:8) = str;
  pstr = ['print -dpng ',filename];
  disp(pstr);
  eval(pstr);
end;

shg;

clear afterframe;
clear mapc2m;
clear mapc2m_inclusions;
