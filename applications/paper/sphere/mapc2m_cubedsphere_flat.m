function [theta,phi,z] = mapc2m_cubedsphere_flat(xc,yc)

[xp,yp,zp] = mapc2m_cubedsphere(xc,yc);

phi = asin(zp);        % in [-pi/2, pi/2]
theta = atan2(yp,xp);  % in [-pi, pi]

mzero = (xp == 0) & (yp == 0);
theta(mzero) = nan;

% Flag = 0 : Return theta in [-pi, pi]
% Flag = 1 : Return theta in [0, pi]
theta = fix_block(theta,1);

%{

tol = pi/2;
switch blockno
    case {0, 3}
        m = abs(theta1 - pi/2) > pi/4 + tol;
        theta1(m) = nan;
    case {6,9}
        % Split 6 over seam at pi
        t = theta1;
        m = t < 0;
        t(m) = t(m) + 2*pi;        
        m = abs(t - pi) > pi/4 + tol;
        theta1(m) = nan;
    case {7,10}
        m = abs(theta1 + pi/2) > pi/4 + tol;
        theta1(m) = nan;        
    case {8,11}
        % Block 0
        m = abs(theta1) > pi/4 + tol;
        theta1(m) = nan;
end
%}
z = 0*theta;

end


function theta = fix_block(theta,flag)

% First, transform theta to [-pi,pi]
t = theta;
m = theta < 0;
t(m) = t(m) + 2*pi;
tmax = max(t(:));
tmin = min(t(:));
if (tmax > pi)
    % Whole block is < 0
    m = theta > 0;
    theta(m) = theta(m)-2*pi;
elseif (tmin < pi)
    % Whole block is > 0
    m = theta < 0;
    theta(m) = theta(m) + 2*pi;
end
tmin = min(theta(:));
tmax = max(theta(:));
if (tmin*tmax < 0)
    % This shouldn't happen, since a block should be
    % entirely in [-pi,0] or [0,pi] after fix
    error('cubedsphere_flag: Theta changes sign');
end
        
if flag == 1
    % return theta in [0,2*pi]
    if (tmax <= 0)
        theta = theta + 2*pi;
    end

end


end