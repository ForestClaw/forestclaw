function [xp,yp,zp] = mapc2m_sphere(xc1,yc1)
%
% Map [0,1] to the hemisphere.  Use the block number to
% determine sign of z-coordinate, i.e. which hemisphere
% should be mapped.
%

% First, get into [-1 1] from [0,1]
xc = 2*xc1 - 1;
yc = 2*yc1 - 1;

r1 = 1;

xp = xc;
yp = yc;
blockno = getblocknumber();

% Map any cells in the ghost cell region to the interior
d = max(xc-1,0) + max(-1-xc,0);
xc = (1-d)./(1+d).*xc;

d = max(yc-1,0) + max(-1-yc,0);
yc = (1-d)./(1+d).*yc;

% Map to the circle
[xp,yp] = mapc2m_disk_sp(xc,yc);

r2 = min(xp.^2 + yp.^2,1);
zp = sqrt(1 - r2);

% Flip the sign on any ghost cell values
mghost = (abs(xc) > 1) | (abs(yc) > 1);

zp(mghost) = -zp(mghost);          % negate z in lower hemisphere

% Set the sign on the grid for the southern hemisphere
blockno = getblocknumber();

if (blockno == 1)
  zp = -zp;
end;

% ----------------------------------------------------------------
% Disk mapping that has been stretched for the sphere
% ----------------------------------------------------------------
function [xp,yp] = mapc2m_disk_sp(xc,yc)

% Disk mapping used for the sphere.

d = max(abs(xc),abs(yc));  % value on diagonal of computational grid
d = max(d, 1e-10);         % to avoid divide by zero at origin

D = sin(pi*d/2)/sqrt(2);   % Smoother map?
R = ones(size(d));         % mapping d to R(d) from (??)

center = D - sqrt(R.^2 - D.^2);

xp = D./d .* abs(xc);
yp = D./d .* abs(yc);

ij = find(abs(yc)>=abs(xc));
yp(ij) = center(ij) + sqrt(R(ij).^2 - xp(ij).^2);

ij = find(abs(xc)>=abs(yc));
xp(ij) = center(ij) + sqrt(R(ij).^2 - yp(ij).^2);

xp = sign(xc) .* xp;
yp = sign(yc) .* yp;
