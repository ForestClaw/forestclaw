function [xp,yp,zp] = mapc2m_latlong2(xc,yc,R,phi0,phi1)

% assume that [xc,yc] in a lat/long box

theta = 2*pi*xc;

% yc in [0,1]
phi = pi*(phi0 + (phi1-phi0)*yc);

xp = R*cos(phi).*cos(theta);
yp = R*cos(phi).*sin(theta);
zp = R*sin(phi);

end
