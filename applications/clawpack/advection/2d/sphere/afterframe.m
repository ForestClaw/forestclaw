setviews;

yrbcolormap;
caxis([0 1]);

showpatchborders;
setpatchborderprops('linewidth',1);
daspect([1 1 1]);
set(gca,'box','on');
view(3);

cv = linspace(0,1,11);
cv([1 end]) = [];
drawcontourlines(cv);

view(vbot)

MaxFrames = 64;
NoQuery = 0;
prt = false;
if (prt)
    filename = sprintf('sphere_%0.4d.png',Frame);
    disp(filename);
    print('-dpng',filename);
end

shg