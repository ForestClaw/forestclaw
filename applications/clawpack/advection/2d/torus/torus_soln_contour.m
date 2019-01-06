function  [xpout,ypout,zpout] = torus_soln_contour(T,alpha,revs_per_sec,period)
%
% TORUS_SOLN_CONTOUR returns location of 0 contour of solution on torus
%
% TORUS_SOLN_CONTOUR(T,alpha,R) returns location of boundary of a 
% disk as it is is advected using a streamfunction
%
%    psi(theta,phi) = (2*pi*R)*alpha*(theta + pi + alpha*sin(theta+phi))
%
% where R is the desired number of revolutions per second and alpha is the
% ratio of the smaller radius to the larger.
%
% The initial "disk" is actually obtained by intersecting a sphere, located
% at (1,0,r0), of r0 is the inner radius of the torus.
% 


global N_torus period_torus alpha_torus revs_per_second_torus;
global r0_torus xc0_torus yc0_torus zc0_torus;

if (nargin < 4)
    period = 0;   % Won't be used.
    if (nargin < 2)
        alpha = 0.4;
        revs_per_sec = 0.5;
        period = 16;
    end
end

alpha_torus = alpha;
revs_per_second_torus = revs_per_sec;
r0_torus = alpha;  % assumes outer radius == 1

period_torus = period;

N = 500;
xc0_torus = 1;    % center of sphere used to initialize problem
yc0_torus = 0;
zc0_torus = alpha_torus;
[x0,y0,z0] = torus_init(N);

N_torus = length(x0);

use_streamfunction = false;

if (T > 0)
    % Solve for t > 0
    N_torus = length(x0);

    if (use_streamfunction)
        % More general than just rigid body rotation
        f0 = [x0; y0; z0];
        opt.RelTol = 1e-12;
        [~,pout] = ode45(@psi_rhs_simple,[0,T],f0,opt);
        
        xp = pout(end,1:N_torus);
        yp = pout(end,(N_torus+1):(2*N_torus));
        zp = pout(end,(2*N_torus+1):end);
    else
        
        % Rigid body rotation (faster than above?)
        
        th = 2*pi*revs_per_sec*T;
        xp = cos(th)*x0 + sin(th)*y0;
        yp = sin(th)*x0 - cos(th)*y0;
        zp = z0;
    end
else
    % Return initial conditions
    xp = x0;
    yp = y0;
    zp = z0;
end

if (nargout > 0)
    % We are plotting the contour on top of a ForestClaw simulation
    xpout = xp;
    ypout = yp;
    zpout = zp;
else
    % This function is being called as a standalone function.
    plot_torus(N);
    hold on
    
    o = findobj('tag','torus_curve');
    if ~isempty(o)
        delete(o)
    end
    
    pline = plot3(xp(:),yp(:),zp(:),'k-','linewidth',2,'markersize',20);
    set(pline,'tag','torus_curve')
end

end

function plot_torus(N)

o = findobj('tag','torus');
if (~isempty(o))
    delete(o);
end

xc = linspace(0,1,N);
yc = linspace(0,1,N);
[xcm,ycm] = meshgrid(xc,yc);

[xp,yp,zp] = mapc2m_torus(xcm,ycm);
p = surf(xp,yp,zp);
set(p,'edgecolor','none');
set(p,'facecolor',0.8*[1,1,1]);
set(p,'facealpha',1);
daspect([1,1,1])
set(p,'tag','torus');

o = findobj('type','light');
if ~isempty(o)
    delete(o);
end
camlight

end


% Get points that make up zero-level contour
function [xdata,ydata,zdata] = torus_init(N)

global alpha_torus

xc = linspace(0,1,N);
yc = linspace(0,1,N);
[xcm,ycm] = meshgrid(xc,yc);

fm = torus_surface(xcm,ycm);

h = contourc(xc,yc,fm,[0,0]);

m = length(h);
xdata = nan(m,1);   % Final length < m
ydata = nan(m,1);
zdata = nan(m,1);

st_idx = 1;
k = 1;
while (1)
  % cval = h(1,st_idx);
  next_length = h(2,st_idx);

  xl = h(1,st_idx+1:st_idx+next_length);
  yl = h(2,st_idx+1:st_idx+next_length);

  [x0,y0,z0] = mapc2m_torus(xl(:),yl(:), alpha_torus);

  xdata(k:k+next_length-1) = x0;
  ydata(k:k+next_length-1) = y0;
  zdata(k:k+next_length-1) = z0;
  k = k + next_length;
  st_idx = st_idx + next_length + 1;
  if (st_idx > length(h))
      break;
  end
end

% Remove extra nan's at the end of the array.
m = isnan(xdata);
xdata(m) = [];
ydata(m) = [];
zdata(m) = [];


end

function f = torus_surface(xc,yc)

global r0_torus xc0_torus yc0_torus zc0_torus alpha_torus;

[xp,yp,zp] = mapc2m_torus(xc,yc,alpha_torus);


r2 = (xp-xc0_torus).^2 + (yp-yc0_torus).^2 + (zp - zc0_torus).^2 - r0_torus^2;
f = r2;

end

% Is this computed from a streamfunction?
function fp = psi_rhs(t,f)

global N_torus period_torus alpha_torus revs_per_second_torus;

R = revs_per_second_torus;
N = N_torus;

x = f(1:N)';
y = f((N+1):(2*N))';
z = f((2*N+1):end)';

r = sqrt(x.^2 + y.^2);

t1 = [-y; x; zeros(size(x))];
t2 = [-z.*x./r; -z.*y./r; r - 1];

sinphi = z/alpha_torus;
cosphi = (r-1)/alpha_torus;

% cos(theta + phi)
costpp = (x./r).*cosphi - (y./r).*sinphi;

u = ((2*pi)*R./r).*(1 + alpha_torus*costpp).*(t1(1,:)-t2(1,:));
v = ((2*pi)*R./r).*(1 + alpha_torus*costpp).*(t1(2,:)-t2(2,:));
w = ((2*pi)*R./r).*(1 + alpha_torus*costpp).*(t1(3,:)-t2(3,:));

fp = -[u(:); v(:); w(:)];

if (period_torus > 0)
    vt = -cos(2*pi*t/period_torus);
    fp = vt*fp;
end


end

% Get the velocity field directly from velocity field for rigid body
% rotation.
function fp = psi_rhs_simple(t,f)

global N_torus period_torus revs_per_second_torus;

R = revs_per_second_torus;
N = N_torus;

x = f(1:N)';
y = f((N+1):(2*N))';
% z = f((2*N+1):end)';

% r = sqrt(x.^2 + y.^2);

% Covariant vectors
t1 = [-y; x; zeros(size(x))];
% t2 = [-z.*x./r; -z.*y./r; r - 1];

t1t = t1';
% Rigid body rotation : u = (2*pi*r)*revs_per_sec*t1./r = 2*pi*R*t1
fp = (2*pi*R)*t1t(:);

if (period_torus > 0)
    vt = -cos(2*pi*t/period_torus);
    fp = vt*fp;
end

end
