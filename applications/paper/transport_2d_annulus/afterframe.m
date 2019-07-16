global square

setviews;

alpha = 0.4;
s = 1e-2;    
alim = [-1,1];
axis([alim alim]);
daspect([1 1 1]);
view(vtop)

showpatchborders(1:10);
setpatchborderprops('linewidth',1)
caxis([0,1])

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

if (mq == 3)
    % Plot the error
    ca = [-max([qmin,qmax]),max([qmin,qmax])];
    colormap(jet);
else    
    % Plot the solution
    yrbcolormap
    ca = [0, 1];
end

% colorbar;
caxis(ca)
hidegridlines
if square
    axis([0,1,0,1]);
end

% showgridlines;

colormap(parula);
cmax = max(abs([qmin,qmax]));
ca = [-cmax,cmax];
caxis(ca);
caxis([-1,1]);
% caxis([-1,1]*1e-14);

if (square)
    daspect([1,5,1]);
end

set(gca,'fontsize',16);

axis([-0.707106781186547   0.707106781186548   0.282842712474619,1]);

plot_path = false;
if (plot_path)
    th = 2*pi*linspace(0.25-1/32,0.25+1/32,200);
    xpath = ravg*cos(th);
    ypath = ravg*sin(th);
    ravg = 0.7;
    hold on;
    plot(xpath,ypath,'k','linewidth',3);
    plot(xpath([1,end]),ypath([1 end]),'r','linewidth',2);
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
