function test_cons(f)

if (nargin == 0)
    f = 0;
end

R = 1;
H = 2*pi*R;
u = [1,1];

exact = true;


minlevel = 1;
maxlevel = 2;

mx0 = 8;
mx = mx0*2^f;

dt0 = 2e-2/2^f;
dt = dt0/2^(minlevel-1);

% ---------------------------------------
e = @(h) R - sqrt(R^2 - (h/2).^2);

% Area of coarse grid cell
Ag = @(h) 0.5*e(h).*(h/2);

% Length of cell on fine grid
len = @(h) sqrt(e(h).^2 + (h/2).^2);

% Curvature in cylinder
kappa = @(h) -(2*len(h) - h)./(2*Ag(h));

% ---------------------------------------

Nc = mx*2^minlevel;
Nf = mx*2^maxlevel;

dthc = 2*pi/Nc;
dthf = 2*pi/Nf;

if (exact)
    dxc = R*dthc;
    dxf = R*dthf;
else
    dxc = 2*R*sin(dthc/2);
    dxf = 2*R*sin(dthf/2);
end

hcz = H/Nc;
hfz = H/Nf;

Af = dxf*hfz;
% Ac = dxc*hcz;
Ac = Af*4^(maxlevel-minlevel);  % Coarse grid area is approx. as 4*(Af)

un = 2*pi*R*u(2);   % f(Q) dot n = Q*(u dot n);  Assumes Q = 1
area_gap = Ag(dxc);
k = kappa(dxc);
delta = -2*area_gap*k*un;    % un = f(Q) dot n_s
cl = dt*delta/Ac;   % mystery factor of 2*pi? 

level = 1;
rf = 2^(maxlevel-level);
N = mx*2^(level)*rf;
% fprintf('%5d %5d %5d %24.16f %24.16f %24.16f\n',minlevel,rf,N,dxf,hfz,Ac);


level = 2;
rf = 2^(maxlevel-level);
N = mx*2^(level)*rf;
% fprintf('%5d %5d %5d %24.16f %24.16f %24.16f\n',maxlevel,rf,N,dxf,hfz,Af);

fprintf('%20s %8d\n','mx',mx);
fprintf('%20s %8d\n','Nf',Nf);
fprintf('%20s %8d\n','Nc',Nc);
fprintf('\n');
fprintf('%20s %24.16f\n','dxc',dxc);
fprintf('%20s %24.16f\n','hcz',hcz);
fprintf('%20s %24.16f\n','dxf',dxf);
fprintf('%20s %24.16f\n','hfz',hfz);
fprintf('\n');
fprintf('%20s %24.16e\n','Cell area (coarse)',Ac);
fprintf('%20s %24.16e\n','Cell area (fine)',Af);
fprintf('\n');
fprintf('%20s %24.16e\n','curvature',k);
fprintf('\n');
fprintf('%20s %24.16e\n','dt*delta',dt*delta);
fprintf('\n');

area_total = (Ac*Nc*Nc + Af*Nf*Nf)/2;

fprintf('%20s %24.16e\n','Total area',area_total);

fprintf('%20s %24.16e\n','Loss of conservation',Nc*cl*Ac);


end