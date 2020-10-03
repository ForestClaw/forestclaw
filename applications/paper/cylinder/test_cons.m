function test_cons(f)

global R H exact phi0 phi1 f_phi

if (nargin == 0)
    f = 0;
end

map = 'latlong';    % 'cylinder', 'latlong'

switch map
    case {'cylinder'}
        R = 1;
        H = 2*pi*R;
        u = [1,1];
        
        % Basis vectors
        tensors = @cylinder_tensors;
        
    case {'latlong'}
        R = sqrt(pi^2 + 1);
        arc = pi/2 - acos(pi/R);
        phi0 = -0.5*arc/(pi/2);
        phi1 = -phi0;
        u = [1,1];
        tensors = @latlong_tensors;
end

exact = false;
f_phi = 0.95;   % Fraction from 0 to 1;  0.5 is the equator

minlevel = 1;
maxlevel = 2;

mx0 = 8;
mx = mx0*2^f;

dt0 = 2e-2/2^f;
dt = dt0/2^(minlevel-1);

%%{
% ---------------------------------------
e = @(h,R) R - sqrt(R^2 - (h/2).^2);

% Area of coarse grid cell
Ag = @(h,R) 0.5*e(h,R).*(h/2);

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
% dxf = dxc/2^(maxlevel-minlevel);
% dyf = dyc/2^(maxlevel-minlevel);


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

H_gap = ec*nc - (ef1*nf1 + ef2*nf2);  % Vector : -2*Ac*k*ns
delta = dot(H_gap,uvel);
cl = dt*delta/Ac;   % mystery factor of 2*pi? 

fprintf('%20s %8d\n','mx',mx);
fprintf('%20s %8d\n','Nf',Nf);
fprintf('%20s %8d\n','Nc',Nc);
fprintf('\n');
ns = H_gap/norm(H_gap,2);
fprintf('%20s (%.8f, %.8f, %.8f)\n','Ns',ns(1),ns(2),ns(3));
fprintf('\n');
phi = pi*(phi0 + (phi1 - phi0)*f_phi);
h = 2*R*cos(phi)*sin(pi*dxc);
area_gap = Ag(h,cos(phi)*R);
k = norm(H_gap,2)/(2*area_gap);
fprintf('%20s %24.16f\n','curvature',k);
fprintf('%20s %24.16f\n','True curvature',1/R);
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
cc = mdpt + dc*mdir/4;
x = cc(1);
y = cc(2);
[t1,t2] = cylinder_basis(x,y);

% Compute normals at coarse grid edge
x = mdpt(1);
y = mdpt(2);
[~,~,n1,~] = cylinder_basis(x,y);
nc = n1;   % Normal to bottom edge

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


function [ep1, ep2, mdir] = latlong_endpoints(dx,dy)

global f_phi;

% 0.5 is the equator
ep1 = [0,f_phi];
ep2 = [dx,f_phi];

mdir = [0,-1];   % Direction of coarse grid center from midpoint of edge

end

function [t1,t2,nc,nf1,nf2,ec,ef1,ef2,Ac,Af] = latlong_tensors(dx,dy)

global R exact phi0 phi1

[ep1,ep2,mdir] = latlong_endpoints(dx,dy);

% Midpoint and length of edge
mdpt = (ep1 + ep2)/2;
d = (ep2-ep1)/norm(ep2-ep1,2);
dc = norm(ep1-ep2,2);

% Coarse grid center (hc = dc/2;  center is half way down
cc = mdpt + dc*mdir/4;
x = cc(1);
y = cc(2);
[t1,t2] = latlong_basis(x,y);

% Compute normals at coarse grid edge
x = mdpt(1);
y = mdpt(2);
[~,~,~,n2] = latlong_basis(x,y);
nc = n2;   % Normal to bottom edge

pf1 = ep1 + (dc/4)*d;
x = pf1(1);
y = pf1(2);
[~,~,~,n2] = latlong_basis(x,y);
nf1 = n2;   % Normal to bottom edge

pf1 = ep1 + (3*dc/4)*d;
x = pf1(1);
y = pf1(2);
[~,~,~,n2] = latlong_basis(x,y);
nf2 = n2;   % Normal to bottom edge

% Get edge lengths
dth = 2*pi*dx;
fe = pf1 + dc/4*d;
fc = fe - (dc/4)*mdir;
ye = fe(2);
phi = pi*(phi0 + (phi1-phi0)*ye);
dphi = pi*(phi1-phi0)*dy;
Rphi = R*cos(phi);
if (exact)
    ec = Rphi*dth;
    ef1 = ec/2;   
    Af = R^2*cos(phi)*(2*pi*dx)*dphi;
else
    ec = 2*Rphi*sin(dth/2);
    ef1 = 2*Rphi*sin(dth/4);
    ez = 2*R*sin(dth/4);
    Af = ef1*ez;
end
ef2 = ef1;

% Areas
Ac = 4*Af;

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
