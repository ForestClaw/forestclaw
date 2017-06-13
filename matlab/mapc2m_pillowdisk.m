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


maptype = 'disk';

% Get coordinates into [-1,1]x[-1,1] 
[xc,yc,zc] = mapc2m_cart(xc1,yc1);


xp = xc;
yp = yc;
zp = 0*xp;

r1 = 1;

% ijlower = find(xc < -1);           % indices of points in lower hemisphere
% xc(ijlower) = -2 - xc(ijlower);

d = max(abs(xc),abs(yc));  % value on diagonal of computational grid
d = max(d, 1e-10);         % to avoid divide by zero at origin

switch maptype
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
  zp = sqrt(r1^2 - (xp.^2 + yp.^2));
  zp(ijlower) = -zp(ijlower);          % negate z in lower hemisphere
end;


if (exist('rrot.dat'))
  rrot = load('rrot.dat');

  [m,n] = size(xp);
  xpv = reshape(xp,1,m*n);
  ypv = reshape(yp,1,m*n);
  zpv = reshape(zp,1,m*n);

  v = rrot*[xpv; ypv; zpv];
  xp = reshape(v(1,:),m,n);
  yp = reshape(v(2,:),m,n);
  zp = reshape(v(3,:),m,n);
end;
