function test_cons(f)

global R H exact phi0 phi1 f_phi exact_discrete

if (nargin == 0)
    f = 0;
end

map = 'cylinder';    % 'cylinder', 'latlong'

switch map
    case {'cylinder'}
        R = 1;
        H = 2*pi*R;
        u = [1,1];
        
        % Basis vectors
        tensors = @cylinder_tensors;
        
    case {'latlong'}
        % Construct sphere by "blowing out" the cylinder, so that 
        % cylinder fits exactly inside.  
        R = sqrt(pi^2 + 1);
        arc = pi/2 - acos(pi/R);
        phi0 = -0.5*arc/(pi/2);
        phi1 = -phi0;
        u = [1,1];
        tensors = @latlong_tensors;
end

exact = false;
exact_discrete = false;  % For spherical coordinates
f_phi = 0.65;   % Fraction from 0 to 1;  0.5 is the equator

minlevel = 1;
maxlevel = 2;

mx0 = 8;
mx = mx0*2^f;

dt0 = 2e-2/2^f;
dt = dt0/2^(minlevel-1);

%%{
% ---------------------------------------
e = @(h,R) R - sqrt(R^2 - (h/2).^2);

% Area of the triangular gap
Ag = @(h,R) e(h,R).*h/2;

% Length of cell on fine grid
len = @(h,R) sqrt(e(h,R).^2 + (h/2).^2);

% Curvature in cylinder
kappa = @(h,R) -(2*len(h,R) - h)./(2*Ag(h,R));

% ---------------------------------------
%}

Nc = mx*2^minlevel;
Nf = mx*2^maxlevel;

dxc = 1/Nc;
dyc = 1/Nc;


% Get tensors at places needed to compute velocity, edge lengths and 
% normals
[t1,t2,nc,nf1,nf2,ec,ef1,ef2,Ac,Af] = tensors(dxc,dyc);


% Velocity in coarse cell center needed to compute f(Q)
% Choose first cell centered point in coarse cell
uvel = u(1)*t1 + u(2)*t2;  

% un = 2*pi*R*u(2);   % f(Q) dot n = Q*(u dot n);  Assumes Q = 1
% area_gap = Ag(dxc);
% k = kappa(dxc);
% delta = -2*area_gap*k*un;    % un = f(Q) dot n_s

H_gap = ec*nc - (ef1*nf1 + ef2*nf2);  % Vector : -2*Ag*k*ns
delta = dot(H_gap,uvel);
cl = dt*delta/Ac;   % mystery factor of 2*pi? 

fprintf('%20s %8d\n','mx',mx);
fprintf('%20s %8d\n','Nf',Nf);
fprintf('%20s %8d\n','Nc',Nc);
fprintf('\n');
ns = H_gap/norm(H_gap,2);
fprintf('%20s (%.8f, %.8f, %.8f)\n','Ns',ns(1),ns(2),ns(3));
fprintf('\n');
if (strcmp(map,'latlong'))
    phi = pi*(phi0 + (phi1 - phi0)*f_phi);
    dth = 2*pi*dxc;
    h = 2*R*cos(phi)*sin(dth/2);
    area_gap = Ag(h,R*cos(phi));
    k = norm(H_gap,2)/(2*area_gap);
    true_kappa = 1/R;
else
    dth = 2*pi*dxc;
    h = 2*R*sin(dth/2);
    area_gap = Ag(h,R);
    k = norm(H_gap,2)/(2*area_gap);
    k = kappa(h,R);
    true_kappa = 0.5/R;
end
fprintf('%20s %24.16f\n','curvature',k);
fprintf('%20s %24.16f\n','True curvature',true_kappa);
fprintf('%20s %24.8e\n','Curvature (error)',abs(true_kappa - k));
fprintf('\n');
fprintf('%20s %24.16f\n','dxc',ec);
fprintf('%20s %24.16f\n','dxf',ef1);
fprintf('\n');
fprintf('%20s %24.16e\n','Cell area (coarse)',Ac);
fprintf('%20s %24.16e\n','Cell area (fine)',Af);
fprintf('\n');
% fprintf('\n');
fprintf('%20s %24.16e\n','dt*delta',dt*delta);
fprintf('\n');

if (strcmpi(mapping,'cylinder') == 1)
    area_total = (Ac*Nc*Nc + Af*Nf*Nf)/2;
    fprintf('%20s %24.16e\n','Total area',area_total);
    fprintf('\n');
end

fprintf('%20s %24.16e\n','Loss of conservation',Nc*cl*Ac);


end

function [ep1, ep2, mdir] = cylinder_endpoints(dx,~)

ep1 = [0,0.5];
ep2 = [dx,0.5];

mdir = [0,-1];   % Direction of coarse grid center from midpiont of edge

end

function [t1,t2,nc,nf1,nf2,ec,ef1,ef2,Ac,Af] = cylinder_tensors(dx,dy)

global H R exact

[ep1,ep2,mdir] = cylinder_endpoints(dx,dy);

% Midpoint and length of edge
mdpt = (ep1 + ep2)/2;
dc = norm(ep1-ep2,2);

% Coarse grid center (hc = dc/2;  center is half way down
cc = mdpt + dc*mdir/2;
x = cc(1);
y = cc(2);
[t1,t2] = cylinder_basis(x,y);

% Compute normals at coarse grid edge
x = mdpt(1);
y = mdpt(2);
[~,~,~,n2] = cylinder_basis(x,y);
nc = n2;   % Normal to bottom edge

x = mdpt(1) + dc/4;
y = mdpt(2);
[~,~,~,n2] = cylinder_basis(x,y);
nf1 = n2;   % Normal to bottom edge

x = mdpt(1) + 3*dc/4;
y = mdpt(2);
[~,~,~,n2] = cylinder_basis(x,y);
nf2 = n2;   % Normal to bottom edge

% Get area of fine grid
dth = 2*pi*dx;
if (exact)
    hxc = R*dth;
    hxf = R*dth/2;
else
    hxc = 2*R*sin(dth/2);
    hxf = 2*R*sin(dth/4);
end
hyf = H*dy/2;
Af = hxf*hyf;
Ac = 4*Af;

ec = hxc;
ef1 = hxf;
ef2 = hxf;

end


function [t1,t2,nx,ny] = cylinder_basis(x,y)

global R H

% Compute tangents and normals at specified (x,y) location

th = 2*pi*x;
t1 = 2*pi*R*[-sin(th); cos(th); 0];
t2 = [0;0;H];

if (nargout > 2)
    G = [dot(t1,t1), dot(t1,t2); dot(t2,t1), dot(t2,t2)];

    Ginv = inv(G);

    t1inv = Ginv(1,1)*t1 + G(1,2)*t2;
    t2inv = Ginv(2,1)*t1 + G(2,2)*t2;

    nx = t1inv./norm(t1inv,2);
    ny = t2inv./norm(t2inv,2);
end

end


function [ep1, ep2, ep3, mdir] = latlong_endpoints(dx,dy)

global f_phi;

% Endpoints are in (theta,phi) coordinates.
% 0.5 is the equator
ep1 = [0,f_phi];
ep2 = [dx,f_phi];
mdpt = (ep1 + ep2)/2;

mdir = [0,-1];   % Direction of coarse grid center from midpoint of edge

% Location of coarse grid center
ep3 = mdpt + dy*mdir/2;

end

function [t1,t2,nc,nf1,nf2,ec,ef1,ef2,Ac,Af] = latlong_tensors(dx,dy)

global R exact phi0 phi1 exact_discrete

mysinc = @(x) sin(x)./x;

% Return endpoints of coarse grid cell, along with
% direction to coarse grid center. 
[ep1,ep2,ep3,mdir] = latlong_endpoints(dx,dy);

% Midpoint and length of edge in (x,y) coordinates
mdpt = (ep1 + ep2)/2;
dxc = norm(ep1-ep2,2);
d1 = (ep2-ep1)/dxc;
% d2 = mdir;

% Compute tangents at center of coarse grid cell (needed for computing
% velocity). 
x = ep3(1);
y = ep3(2);
[t1,t2] = latlong_basis(x,y);

% Assume phi is constant along coarse edge
phi = pi*(phi0 + (phi1-phi0)*mdpt(2));
dphi = pi*(phi1-phi0)*dy;

theta = 2*pi*mdpt(1);  % x in [0,1]
dth = 2*pi*dx;  % equal to dc ? 

% pf1 = ep1 + (3*hc/4)*d;
% fe = pf1 + dc/4*d;
% %fc = fe - (dc/4)*mdir;  % center of fine grid ghost cell
% ye = fe(2);
Rphi = R*cos(phi);

if (exact)
    Rphi = R*cos(phi);
    
    % Coarse grid edge lengths and normal
    nc_average = [-sin(phi)*cos(theta)*mysinc(dth/2)
                  -sin(phi)*sin(theta)*mysinc(dth/2)
                  cos(phi)];
    lcnc = Rphi*dth*nc_average;
    ec = norm(lcnc,2);
    nc = lcnc/ec;
    
    % Compute fine grid edge lengths and normal
    b = (1-cos(dth/2))/(dth/2);
    w = [sin(phi)*sin(theta)
         -sin(phi)*cos(theta)
         0];
        
    nf1_average = nc_average - b*w;
    lfnf1 = Rphi*(dth/2)*nf1_average;
    ef1 = norm(lfnf1,2);
    nf1 = lfnf1/ef1;  

    nf2_average = nc_average + b*w;
    lfnf2 = Rphi*(dth/2)*nf2_average;
    ef2 = norm(lfnf2,2);
    nf2 = lfnf2/ef2;  

    if (abs(ec*nc - (ef1*nf1 + ef2*nf2)) > 1e-12) 
        error('Something went wrong with exact latlong formula');
    end
    
    % Exact area of fine grid
    Af = R^2*cos(phi)*dth*dphi;
    Ac = 4*Af;
else
    % Compute discrete normals and edge lengths
    
    % Fine grid edge normals
    pf1 = ep1 + (dxc/4)*d1;
    x = pf1(1);
    y = pf1(2);
    [~,~,~,n2] = latlong_basis(x,y);
    nf1 = n2;   % Normal to bottom edge

    pf1 = ep1 + (3*dxc/4)*d1;
    x = pf1(1);
    y = pf1(2);
    [~,~,~,n2] = latlong_basis(x,y);
    nf2 = n2;   % Normal to bottom edge
    
    % Get fine grid edge lengths
    ef1 = 2*Rphi*sin(dth/4);
    ef2 = 2*Rphi*sin(dth/4);
    
    % Compute coarse grid
    if (exact_discrete)
        lcnc = ef1*nf1 + ef2*nf2;
        ec = norm(lcnc,2);
        nc = lcnc/ec;
    else
        % Compute normals at coarse grid edge
        x = mdpt(1);
        y = mdpt(2);
        [~,~,~,n2] = latlong_basis(x,y);
        nc = n2;   % Normal to bottom edge    
        ec = 2*Rphi*sin(dth/2);
    end    
    ez = 2*R*sin(dth/4);
    Af = ef1*ez;
    
    % Areas
    Ac = 4*Af;
end

end

function [t1,t2,nx,ny] = latlong_basis(x,y)

global R phi0 phi1

% Compute tangents and normals at specified (x,y) location

th = 2*pi*x;
dthdx = 2*pi;
phi = pi*(phi0 + (phi1 - phi0)*y);
dphidy = pi*(phi1 - phi0);
t1 = R*cos(phi)*[-sin(th); cos(th); 0]*dthdx;
t2 = R*[-sin(phi)*cos(th); -sin(phi)*sin(th); cos(phi)]*dphidy;

if (nargout > 2)
    G = [dot(t1,t1), dot(t1,t2); dot(t2,t1), dot(t2,t2)];

    Ginv = inv(G);

    t1inv = Ginv(1,1)*t1 + G(1,2)*t2;
    t2inv = Ginv(2,1)*t1 + G(2,2)*t2;

    nx = t1inv./norm(t1inv,2);
    ny = t2inv./norm(t2inv,2);
end

end
