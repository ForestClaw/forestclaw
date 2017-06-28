axis([0 1 0 1]);
rybcolormap
caxis([-.1 .1])
colorbar

showpatchborders
setpatchborderprops('linewidth',1);

hold on
plot([0.5 0.5],[0 1],'w','linewidth',3);
hold off
zl = zlim;
set(gca,'zlim',[zl(1) 0]);

axis square

clear afterframe;
