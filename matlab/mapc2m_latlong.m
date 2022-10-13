function [xp,yp,zp] = mapc2m_latlong(xc,yc)

% assume that [xc,yc] in a lat/long box

xp = cosd(yc).*cosd(xc);
yp = cosd(yc).*sind(xc);
zp = sind(yc);

end


