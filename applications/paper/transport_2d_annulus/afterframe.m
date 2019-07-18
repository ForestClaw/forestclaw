global square

setviews;

beta = 0.4;
s = 1e-2;    
alim = [-1,1];
axis([alim alim]);
daspect([1 1 1]);
view(vtop)

showpatchborders(1:10);
setpatchborderprops('linewidth',1)


caxis([-1,1]);

set(gca,'fontsize',16);

axis([-0.707106781186547   0.707106781186548   0.282842712474619,1]);

showgridlines;

plot_path = true;
if (plot_path)
    th = 2*pi*linspace(0.25-1/32,0.25+1/32,200);
    ravg = (1+beta)/2;
    xpath = ravg*cos(th);
    ypath = ravg*sin(th);
    hold on;
    plot(xpath,ypath,'k','linewidth',3);
    plot(xpath([1,end]),ypath([1 end]),'k','linewidth',2);
    hold off;
end


%
NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 8;
  axis([0 1 0 1]);
  filename = sprintf('annulus_%04d.png',Frame)
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
