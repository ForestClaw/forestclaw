global square;

setviews;

if (PlotType == 4)
    axis([0,sqrt(2),-2,2]);
    set(gca,'box','on');
    ca = 1.1*[-max([qmin,qmax]),max([qmin,qmax])];
    % ylim(ca);    
    fprintf('qmin = %24.16e\n',qmin);
    fprintf('qmax = %24.16e\n',qmax);
    shg
    return
end


alpha = 0.4;
s = 1e-2;    
alim = [-1,1];
axis([-1,1 0,1]);
daspect([1 1 1]);
view(vtop)

% yrbcolormap;

showpatchborders(1:10);
setpatchborderprops('linewidth',1)

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

if (mq == 3)
    % Plot the error
    cmax = max(abs([qmin,qmax]));
    ca = [-cmax,cmax];
end

colormap(parula);
cmax = max(abs([qmin,qmax]));
ca = [-cmax,cmax];
caxis(ca);
caxis([-1,1])
% colorbar
% caxis([-1,1]*1e-12);

axis([-0.707106781186547   0.707106781186548   0.282842712474619,1]);

th = 2*pi*linspace(0.25-1/32,0.25+1/32,200);
xpath = ravg*cos(th);
ypath = ravg*sin(th);
ravg = 0.7;
hold on;
plot(xpath,ypath,'k','linewidth',2);
plot(xpath([1,end]),ypath([1 end]),'k','linewidth',2);
hold off;


% showgridlines;
showpatchborders;

if square
    axis([0,1,0,1]);
end

% hidegridlines
% hidepatchborders

% figure(2)
% plot(ycenter,q(1,:),'.-','markersize',20);
% ylim([-1,2]);
% shg
% input('Hit enter to continue : ');

figure(1)


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
