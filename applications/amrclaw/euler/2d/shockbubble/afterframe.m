colormap(jet)

axis([0 2.0 0 0.5])
axis equal; 
axis tight

fprintf('%10s : %12.4f\n','qmin',qmin);
fprintf('%10s : %12.4f\n','qmax',qmax);

caxis([0.1 2.81]);

showpatchborders
setpatchborderprops('linewidth',1);
showgridlines(1:5);
% hidepatchborders;
delete(get(gca,'title'))
set(gca,'fontsize',16);
hidegridlines;

prt = false;
NoQuery = 0;
if (prt)
    MaxFrames = 31;
    fname = sprintf('sbamrmesh%2.2d.png',Frame);
    disp(fname);
    print('-dpng',fname);
end

clear afterframe
clear mapc2m
