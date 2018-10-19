function map_debug2(blockno,center)

figure(3)

if (nargin < 2)
    center = load('center.dat');
end


mx = 8;
h = 1/mx;
mbc = 3;

set_blocknumber(blockno);

xe = linspace(-mbc*h,1+mbc*h,mx+2*mbc+1);
ye = xe;

[xem,yem] = meshgrid(xe,ye);

c = {'r','b','k','m'};

[xepm,yepm,zepm] = mapc2m_bilinear(xem,yem,center);


p = patch(surf2patch(xepm,yepm,zepm));
set(p,'facecolor','none');
set(p,'edgecolor',c{blockno+1});
daspect([1,1,1]);

hold on;
xc = xe(1:end-1) + h/2;
yc = xc;

[xcm,ycm] = meshgrid(xc,yc);

[xcpm,ycpm,zcpm] = mapc2m_bilinear(xcm,ycm,center);

plot(xcpm,ycpm,[c{blockno+1},'o'],'markersize',2);

[xec,yec] = mapc2m_bilinear([0,1,1,0,0],[0,0,1,1,0],center);
plot(xec,yec,'r','linewidth',2);

s = 5*h;
axis([-1-s,1+s,-1-s,1+s]);

shg



end