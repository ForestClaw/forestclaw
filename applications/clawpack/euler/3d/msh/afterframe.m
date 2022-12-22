global region_msh

if (t > 0)
  clim([qmin,qmax]);
end
% cv = linspace(qmin,qmax,30);
% drawcontourlines(cv);
daspect([1 1 1]);
setslicecolor('none','z');

plot_msh();

showslices()
% hideslices('z')

clim([1e3,1e6])
showgridlines
set(gca,'zlimmode','auto')
hidegridlines(5)

view(-17.91,14.98)
set(gca,'clipping','off')


colormap(parula)
colorbar




% hideslices x;
% hideslices y;