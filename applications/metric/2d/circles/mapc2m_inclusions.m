function [xp,yp] = mapc2m_inclusions(xc,yc)

global boxes;

mx1 = xc < -1;
my1 = yc < -1;

mx2 = xc > 1;
my2 = yc > 1;

xc0 = xc;
yc0 = yc;

% isperiodic = false;
% if (isperiodic)
%   xc0(mx1) = xc(mx1) + 2;
%   xc0(mx2) = xc(mx2) - 2;
% 
%   yc0(my1) = yc(my1) + 2;
%   yc0(my2) = yc(my2) - 2;
% end;

xp = xc0;
yp = yc0;

for i = 1:length(boxes),
  % size of grid cells on 'layout' grid

  box = boxes(i);

  xlow = box.xlow;
  ylow = box.ylow;
  xhi = box.xhi;
  yhi = box.yhi;
  r1 = box.r1;
  n = box.size;

  xc1 = -1 + 2*(xc0 - xlow)/(xhi-xlow);
  yc1 = -1 + 2*(yc0 - ylow)/(yhi-ylow);

  square1 = abs(xc1) <= 1 & abs(yc1) <= 1;

  [xp1,yp1] = mapc2m_inclusion(xc1(square1),yc1(square1),r1);

  xp(square1) = (xp1+1)*(xhi-xlow)/2 + xlow;
  yp(square1) = (yp1+1)*(yhi-ylow)/2 + ylow;
end;

% if (isperiodic)
%   xp(mx1) = xp(mx1) - 2;
%   xp(mx2) = xp(mx2) + 2;
%   yp(my1) = yp(my1) - 2;
%   yp(my2) = yp(my2) + 2;
% end;

zp = 0*xp;


function [xp,yp,zp] = mapc2m_inclusion(xc,yc,r1)

%
% MAPC2M_INCLUSION maps the grid [-1,1]x[-1,1] to a grid enclosing a
% circle of radius r1 < 1.
%

xp = xc;
yp = yc;
zp = 0*xc;

d = max(abs(xc),abs(yc));

% -------------------------------------------
% Mapping for inside the circle
% -------------------------------------------
incircle = d <= r1;  % portion of grid deformed for this circle
xc1 = xc(incircle);
yc1 = yc(incircle);
[xp(incircle),yp(incircle)] = mapc2m_incircle(xc1,yc1,r1);


% -------------------------------------------
% Mapping for outside the circle
% -------------------------------------------
outcircle = ~incircle;
xc1 = xc(outcircle);
yc1 = yc(outcircle);
[xp(outcircle),yp(outcircle)] = mapc2m_outcircle(xc1,yc1,r1);


% ------------------------------------------------------------
% Mapping for inside the circle
% ------------------------------------------------------------
function [xp,yp] = mapc2m_incircle(xc,yc,r1)

d = max(abs(xc),abs(yc));
r = max(sqrt(xc.^2 + yc.^2),1e-10);

scale = d./r;

xp = scale.*xc;
yp = scale.*yc;

use_convex = true;
if (use_convex)
  w = (d/r1).^2;
  xp = w.*xp + (1-w).*xc/sqrt(2);
  yp = w.*yp + (1-w).*yc/sqrt(2);
end;

% ---------------------------------------------------------------
% Mapping outside of circle
% ---------------------------------------------------------------
function [xp,yp] = mapc2m_outcircle(xc,yc,r1);

d = max(abs(xc),abs(yc));
r = max(sqrt(xc.^2 + yc.^2),1e-10);

dinv = 1./d;
ct = d./r;

s = (1 - r1*ct)/(1 - r1);
scale = s.*(1 - dinv) + dinv;

xp = scale.*xc;
yp = scale.*yc;
