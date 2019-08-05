function plot_wave_path()

close all;

% Geometric terms
vcart = 1.092505803290319;
tfinal = 0.25;
dlen = tfinal*vcart;
A = 0.01;
rinit = 0.05;
m = 0.5;

beta = 0.4;

ravg = (1 + beta)/2;
thc = pi/2*(1 + 1/8);
pstart = [ravg*cos(thc),ravg*sin(thc)];

% Plot annulus
mi = 4;
mj = 2;
mx = 8;
xc = linspace(0,1,mi*mx + 1);
yc = linspace(0,1,mj*mx + 1);
[xm,ym] = meshgrid(xc,yc);
[xp,yp,zp] = mapc2m(xm,ym);

p = patch(surf2patch(xp,yp,zp));
set(p,'edgecolor','k');
set(p,'facecolor','none');

hold on;

% Plot path
t = linspace(0,tfinal,200);
s = t/tfinal;

xpath = pstart(1) + s*dlen;
ypath = pstart(2) + A*sin(2*pi*m*s);

p_path = plot(xpath,ypath,'r','linewidth',2);
hold on;

% Plot endpoints
pend = pstart + [dlen, 0];
plot(pstart(1),pstart(2),'k.','markersize',30);
plot(pend(1),pend(2),'k.','markersize',30);
plot([pstart(1),pend(1)],[pstart(2),pend(2)],'k','linewidth',2);


% Figure settings
daspect([1,1,1]);
axis([-0.707106781186547, 0.707106781186548, ...
    0.282842712474619,1]);


% Plot and move circle along path
N = 128;
th = linspace(0,2*pi,N+1);

x0 = rinit*cos(th) + pstart(1);
y0 = rinit*sin(th) + pstart(2);

p = plot(x0,y0,'k','linewidth',2);


M = 50;
t = linspace(0,tfinal,M+1);
s = t/tfinal;
for i = 1:M+1
    xth = x0 + s(i)*dlen;
    yth = y0 + A*sin(2*pi*m*s(i));
    
    set(p,'xdata',xth,'ydata',yth);
    pause(0.1);
end
shg;

end