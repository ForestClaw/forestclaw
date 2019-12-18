function [theta,phi,z] = mapc2m_cubedsphere_flat(xc,yc)

[xp,yp,zp] = mapc2m_cubedsphere(xc,yc);

phi = asin(zp);        % in [-pi/2, pi/2]
theta = atan2(yp,xp);  % in [-pi, pi]

tol = 1e-12;
mzero = (abs(xp) < tol) & (abs(yp) < tol);
theta(mzero) = nan;

% Flag = 0 : Return theta in [-pi, pi]
% Flag = 1 : Return theta in [0, 2*pi]
flag = 1;
theta = fix_block(theta,flag);

z = 0*theta;

end


function theta = fix_block(theta,flag)

% atan2(0,-1) maps to pi; 
% atan2(-realmin,-1) maps to -pi

t = theta;
if (max(t(:)) - min(t(:)) > pi)
    % We have wrapped the theta=pi line.  Assume that most of the block is
    % in [-pi,0], so we subtract 2*pi from the mesh points mapped to +pi,
    % since probably only those points on the line y=0, x < 0 have been
    % mapped to +pi.  
    m = theta >= 0;
    t(m) = t(m) - 2*pi;
end
% Mesh points are all in [-pi,0] or in [0,pi]
if flag == 1
    % return theta in [0,2*pi]
    theta = t + pi;
else
    theta = t;
end


end