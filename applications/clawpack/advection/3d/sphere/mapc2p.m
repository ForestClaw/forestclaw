
function [xp,yp,zp] = mapc2p(xc1,yc1,zc1)

global notpillowsphere;

parms = read_vars();

if parms.example == 1
    map = 'cubedsphere';
elseif parms.example == 2
    map = 'pillowsphere';
end

<<<<<<< HEAD
maxelev = parms.maxelev/2;
=======
% b = getblocknumber();
>>>>>>> develop

switch map
    case 'pillowsphere'
        notpillowsphere = false;
        [xp1,yp1,zp1] = mapc2m_pillowsphere(xc1,yc1);

    case 'cubedsphere'
        [xp1,yp1,zp1] = mapc2m_cubedsphere(xc1,yc1);

end
phi = asin(zp1);          % returns value in [-pi/2, pi/2]
theta = atan2(yp1,xp1);    % returns value in [-pi, pi]
m = theta < 0;
theta(m) = theta(m) + 2*pi;

% Assume zc1 in [0,1]
<<<<<<< HEAD
R = maxelev*zc1 + 1;
=======
R = parms.maxelev*zc1 + 1;
>>>>>>> develop
xp = R.*cos(phi).*cos(theta);
yp = R.*cos(phi).*sin(theta);
zp = R.*sin(phi);

end