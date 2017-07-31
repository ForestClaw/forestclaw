s = 1e-2;
axis([-s 1+s -s 1+s])
daspect([1 1 1]);
axis off;

if (~PlotParallelPartitions)        
    colormap(jet);
else
    if (Frame == 0)
        cm = ppcolors(7);
    end
    colormap(cm);
    ppcolors_colorbar(7);

end

cv = linspace(qmin,qmax,21);
cv([1 end]) = [];
drawcontourlines(cv);
setcontourlineprops('linewidth',1);

setpatchborderprops('linewidth',1);
showpatchborders;


view(2);

NoQuery = 0;
prt = false;
if (prt)
  filename = sprintf('quadrants%03d.png',Frame);
  print('-dpng',filename);
end;

shg

clear afterframe;
