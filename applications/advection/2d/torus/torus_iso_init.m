function [th,phi] = torus_iso_init()
close all;

N = 500;
xc = linspace(0,1,N);
yc = linspace(0,1,N);
[xcm,ycm] = meshgrid(xc,yc);

[xp,yp,zp] = mapc2m_torus(xcm,ycm);
p = surf(xp,yp,zp);
set(p,'edgecolor','none');
daspect([1,1,1])
hold on;

fm = iso(xcm,ycm);

h = contourc(xc,yc,fm,[0,0]);
xcont = h(1,:);
ycont = h(2,:);

l1 = h(2,1);
l2 = h(2,l1+2);
h(:,1) = [];
h(:,l1+2) = [];
xinit = [h(1,1:l1), fliplr(h(1,(l1+2):end))];
yinit = [h(2,1:l1), fliplr(h(2,(l1+2):end))];
[xp,yp,zp] = mapc2m_torus(xinit,yinit);
plot3(xp,yp,zp,'k','linewidth',2);
daspect([1,1,1])

% hold on;
% plot3(xp(1),yp(1),zp(1),'b.','markersize',30);
% 
% m2 = xcont < 0.5;
% xinit = xcont(m2);
% yinit = ycont(m2);
% [xp,yp,zp] = mapc2m_torus(xinit,yinit);
% plot3(xp(3:end),yp(3:end),zp(3:end),'g.','linewidth',2);
% hold on;
% plot3(xp(3),yp(3),zp(3),'c.','markersize',30);


if (nargout > 0)
    th = 0;
    phi = 0;
end

end

function f = iso(xc,yc)

[xp,yp,zp] = mapc2m_torus(xc,yc);

r2 = (xp-1).^2 + yp.^2 + (zp - 0.4).^2 - 0.25;

f = r2;



end