function plot_gap()


set(gca,'DataAspectRatiomode','auto')

ax =[-0.009612881581514,   0.010156179730254,  ...
      0.699945546388827,   0.700001822938437];
  
dx = 0.03;
dy = 0.0003;

xc = 0;

xl = xc - dx/2;
xu = xc + dx/2;
yl = 0.7 - dy/1.7;
yu = 0.7 + dy/12;
dp = [0.02, 0.0002];
ax = [xl, xu, yl, yu];

axis(ax);

showgridlines;
setslicecolor([1,1,1]*0.8);
daspect([40,1,1]);

xlabel('x','fontsize',10);
ylabel('y','fontsize',10);
set(gca,'fontsize',10);

delete(get(gca,'title'));
shg


print('-dpng','gap.png');



end