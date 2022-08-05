setviews;

yrbcolormap;
caxis([0 1]);

showpatchborders;
setpatchborderprops('linewidth',1);
daspect([1 1 1]);
set(gca,'box','on');
view(3);

view(vfront)

MaxFrames = 64;
NoQuery = 0;
prt = false;
if (prt)
    filename = sprintf('sphere_%0.4d.png',Frame);
    disp(filename);
    print('-dpng',filename);
end

shg


clear afterframe
clear mapc2m
