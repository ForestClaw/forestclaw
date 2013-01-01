function [xp,yp,zp] = mapc2m(xc1,yc1)
%
% Map [-3,1] x [-1,1] to surface of a sphere of radius r1
%     [-1,1]x[-1,1]  is mapped to the upper hemisphere.
%     [-3,-1]x[-1,1] is mapped to the lower hemisphere.
% Mapping for Figure 5.1(a) as described in Section 4 of
%    Logically Rectangular Grids and Finite Volume Methods for PDEs
%         in Circular and Spherical Domains,
%    by Donna A. Calhoun, Christiane Helzel, and Randall J. LeVeque,
%    http://www.amath.washington.edu/~rjl/pubs/circles


[xp,yp,zp] = mapc2m_pillow(xc1,yc1);

return;

blockno = getblocknumber();

maptype = 'sphere';

xc = 2*xc1 - 1;
yc = 2*yc1 - 1;

xp = xc;
yp = yc;
zp = 0*xp;

if (strcmp(maptype,'cart') == 1)
  return;
end;

% plot_latlong = true;

r1 = 1;

ijlower = find(xc < -1);           % indices of points in lower hemisphere
xc(ijlower) = -2 - xc(ijlower);

d = max(abs(xc),abs(yc));  % value on diagonal of computational grid
d = max(d, 1e-10);         % to avoid divide by zero at origin

switch maptype
case 'cart'
case 'disk'
  D = r1*d/sqrt(2);            % for disk
case {'hemisphere','sphere'}
  % D = r1 * d.*(2-d)/sqrt(2);   % For sphere
  D = r1 * sin(pi*d/2)/sqrt(2);  % Smoother map?
end;
R = r1 * ones(size(d));    % mapping d to R(d) from (??)

center = D - sqrt(R.^2 - D.^2);

xp = D./d .* abs(xc);
yp = D./d .* abs(yc);

ij = find(abs(yc)>=abs(xc));
yp(ij) = center(ij) + sqrt(R(ij).^2 - xp(ij).^2);

ij = find(abs(xc)>=abs(yc));
xp(ij) = center(ij) + sqrt(R(ij).^2 - yp(ij).^2);

xp = sign(xc) .* xp;
yp = sign(yc) .* yp;

switch maptype
case 'disk'
  zp = 0*xp;
case {'hemisphere','sphere'}
  zp = sqrt(1 - (xp.^2 + yp.^2));
  zp(ijlower) = -zp(ijlower);          % negate z in lower hemisphere
end;

if (blockno == 1)
  zp = -zp;
end;

return;

R = load('rrot.dat');

[m,n] = size(xp);

if (m > n)
  z = [xp  yp  zp];

  Rz = R*z';
  xp = Rz(1,:)';
  yp = Rz(2,:)';
  zp = Rz(3,:)';
else
  z = [xp' yp' zp'];
  Rz = R*z';
  xp = Rz(1,:)';
  yp = Rz(2,:)';
  zp = Rz(3,:)';
end;
