function [xp,yp,zp] = mapc2m_cyl_sphere(xc,yc,xbar,R_cyl,H_cyl)

if (ischar(xbar))
    if strcmpi('sphere',xbar) == 1
        xbar = R_cyl*sqrt(1 + 0.25*(H_cyl/R_cyl)^2);
    elseif strcmpi('cylinder',xbar) == 1
        xbar =  R_cyl;
    end
end

assert(R_cyl >  0, 'mapc2m_cyl2pshere : Rcyl == 0');

beta  = H_cyl/R_cyl;
if (xbar == R_cyl)
    % This will return a perfect cylinder.
    c = beta^2/4;
    epc = my_machine_eps(c);
    
    % Compute eps(c) (in case we don't have function `eps(x)` available)
    fprintf('%10s %24.16e\n','eps(c)',eps(c));
    fprintf('%10s %24.16e\n','epc',epc);
   
    e = sqrt(epc/10);
    xbar = R_cyl*(1 + e);
end

% Start of mapping
beta  = H_cyl/R_cyl;
alpha = xbar/R_cyl;
c = beta^2/4;

% Radius of the circle
R_sphere = -R_cyl/2*(1 - alpha + c/(1-alpha));

% Center of circle : (0, -L).  For xbar > R_cyl, L > 0
L = -R_cyl/2*(1 + alpha + c/(1-alpha));

% Get range for phi : 
if (-L <= R_cyl)
    phi0 = asin(-H_cyl/(2*R_sphere));
    phi1 = -phi0;
else
    phi1 = acos((L+R_cyl)/R_sphere);
    phi0 = -phi1;
end
th = 2*pi*xc;
phi = phi0 + yc*(phi1-phi0);


% Radius and height for manifold for circle centered at (0,-L)
R = R_sphere*cos(phi) - L;    
H = R_sphere*sin(phi);

% Spherical coordinates with modified radius.
xp = R.*cos(th);
yp = R.*sin(th);
zp = H;


end

function epx = my_machine_eps(x)

xa = abs(x);

if (xa > 1)
    for k = 1:300
        xa = xa/2;
        if (xa < 1)
            kmax = k-1;
            break;
        end
    end
    epx = 2^(-52+kmax);
else
    for k = 1:300
        xa = 2*xa;
        if (xa > 1)
            kmax = k;
            break;
        end
    end
    epx = 2^(-52-kmax);
end

end