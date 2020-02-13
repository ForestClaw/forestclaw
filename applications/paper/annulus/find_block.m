function find_block()

% close all;

d = load('mapping.dat');

beta = d(3);
t1 = d(4);
t2 = d(5);
theta = [t1, t2];

x = linspace(0,1,5);
y = linspace(0,1,3);
[xm,ym] = meshgrid(x,y);

[xp,yp,zp] = mapc2m_annulus(xm,ym,beta,theta);

p = patch(surf2patch(xp,yp,zp));
set(p,'edgecolor','k');
set(p,'facecolor','none');

daspect([1,1,1]);
axis([-0.707106781186547   0.707106781186548   0.282842712474619,1]);

hold on;

xlow = 0.5;
ylow = 0.5;

brick_data = load('brick.dat');
mi = brick_data(1,1);
mj = brick_data(1,2);

for b = 0:(mi*mj-1)
    set_blocknumber(b);
    s = 0;
    [xc1,yc1,~] = mapc2m_brick(xlow,ylow,s);
    [xp,yp]= mapc2m_annulus(xc1,yc1,beta,theta);
    str = sprintf('%2d',b);
    text(xp,yp,str,'fontsize',16);
    %pause;
end

end