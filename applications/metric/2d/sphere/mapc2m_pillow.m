function [xp,yp,zp] = mapc2m_pillow(xc1,yc1)
%
% Map [-1,1] to the hemisphere.  Use the block number to
% determine sign of z-coordinate.
%

comp_grid = false;

% First, get into [-1 1] from [0,1]
xc = 2*xc1 - 1;
yc = 2*yc1 - 1;

r1 = 1;

xp = xc;
yp = yc;
blockno = getblocknumber();


if (blockno == 0)
  zp = 0*xp + 1;
else
  zp = 0*xp - 1;
end
if (comp_grid)
    return;
end

% This is probably not necessary since we don't expect to have any ghost cells.
d = max(xc-1,0) + max(-1-xc,0);
xc = (1-d)./(1+d).*xc;

d = max(yc-1,0) + max(-1-yc,0);
yc = (1-d)./(1+d).*yc;

[xp,yp] = mapc2m_disk_sp(xc,yc);

r2 = min(xp.^2 + yp.^2,1);
zp = sqrt(1 - r2);

mghost = (abs(xc) > 1) | (abs(yc) > 1);

zp(mghost) = -zp(mghost);          % negate z in lower hemisphere

if (blockno == 1)
    zp = -zp;
end

% [xp,yp,zp] = rotate_map(xp,yp,zp);


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


% --------------------------------------------
% Rotated mapping
% --------------------------------------------
function [xp,yp,zp] = rotate_map(xp,yp,zp)
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
