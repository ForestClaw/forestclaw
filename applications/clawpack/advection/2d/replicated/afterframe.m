yrbcolormap;
setviews;

global map isflat;

brick = load('brick.dat');
% axis([0 brick(1,1) 0 brick(1,2)])
axis([0 2 0 2]);
caxis([0 1]);

showpatchborders;
setpatchborderprops('linewidth',1)

view(2);
delete(get(gca,'title'));
set(gca,'fontsize',16,'box','on');
axis square

shg

clear afterframe;
clear mapc2m;
