function [r,q] = map1d(xgrid,ygrid,qgrid)

% [~,my] = size(xgrid);

% r = xgrid(round(my/2),:);

% q = qgrid(round(my/2),:);

maxlevel = 6;
mx = 32;

N = mx*2^maxlevel;
h = 1/N;

xe = linspace(0,1,N+1);
xc = xe(1:end-1) + h/2;
yc = 0.5 + 2*h + 0*xc;

q = interp2(xgrid,ygrid,qgrid,xc,yc,'linear',nan);
r = interp2(xgrid,ygrid,xgrid,xc,yc,'linear',nan);

r = r(:);
q = q(:);


end