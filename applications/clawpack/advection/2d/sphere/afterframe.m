setviews;
daspect([1 1 1]);
axis off;

yrbcolormap;
caxis([0 1]);

showpatchborders;
setpatchborderprops('linewidth',1);
hidepatchborders(7);
view(3);

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
