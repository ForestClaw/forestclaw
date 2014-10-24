function test_squareddisk(map_in,blockno)

global map_choice;

if (nargin < 2)
    blockno = 4;
    if (nargin < 1)
        map_in  = 'squareddisk';
    end
end

map_choice = map_in;  % 'squareddisk' or 'pillowfivepatch'

close all;

set_blocknumber(blockno);

mx = 8;
mbc = 3;

minlevel = 2;
maxlevel = 2;

plot_coarse = true;
plot_fine = true;


% Plot fine grid
if (plot_fine)
    mxf = mx*2^maxlevel;
    mbcf = mbc*2^(maxlevel-minlevel);
    dxf = 1/mxf;
    x = linspace(-mbcf*dxf,1+mbcf*dxf,mxf+2*mbcf+1);
    y = x;
    [xm,ym] = meshgrid(x,y);

    [xp,yp,zp] = mapc2m(xm,ym);

    q = imag(xp).^2 + imag(yp).^2;
    fprintf('Imaginary part : %12.4e\n',max(q(:)));
    ph = patch(surf2patch(real(xp),real(yp),zp,q));
    set(ph,'edgecolor','b','facecolor','none');
end


% Plot coarse grid
if (plot_coarse)
    mxc = mx*2^minlevel;
    mbcc = mbc;
    dxc = 1/mxc;
    x = linspace(-mbcc*dxc,1+mbcc*dxc,mxc+2*mbcc+1);
    y = x;
    [xm,ym] = meshgrid(x,y);

    [xp,yp,zp] = mapc2m(xm,ym);
    
    
    q = imag(xp).^2 + imag(yp).^2;
    fprintf('Imaginary part : %12.4e\n',max(q(:)));
    ph = patch(surf2patch(real(xp),real(yp),zp,q));
    set(ph,'edgecolor','r');
    set(ph,'facecolor','none');
end

qmin = min(q(:));
qmax = max(q(:));

set(gcf,'renderer','zbuffer')

caxis([0 2*qmax]);
daspect([1 1 1]);
hold on;
c = linspace(0,1,mxf+1);
xc = [0 c 1 fliplr(c) 0 0];
yc = [0 zeros(size(c)) 1 ones(size(c)) 1 0];
[xp,yp,zp] = mapc2m(xc,yc);
plot(xp,yp,'k','linewidth',2);
[xp,yp,zp] = mapc2m(0,0);
p1 = plot(xp,yp,'g.','markersize',30);
[xp,yp,zp] = mapc2m(1,0);
p2 = plot(xp,yp,'y.','markersize',30);
legend([p1 p2],{'(0,0)','(1,0)'});

s = sprintf('mx = %d; minlevel = %d\n',mx,minlevel);
title(s,'fontsize',18);


shg;

end


        