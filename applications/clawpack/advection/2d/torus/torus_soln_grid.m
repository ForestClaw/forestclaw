function  qtrue = torus_soln_grid(xcenter,ycenter,xedge,yedge,...
    T,alpha,revs_per_sec,period)
%
% TORUS_SOLN_GRID returns true solution on torus
%
% TORUS_SOLN_GRID(T,alpha,R) returns location of boundary of a 
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


global period_torus alpha_torus revs_per_second_torus;
global r0_torus xc0_torus yc0_torus zc0_torus

global g_mapc2m g_psi_rhs

global reverse_torus

if (nargin < 4)
    period = 0;   % Won't be used.
    if (nargin < 2)
        alpha = 0.4;
        revs_per_sec = 0.5;
        period = 16;
    end
end

g_psi_rhs = @psi_rhs_torus;
g_mapc2m = @mapc2m_torus;

alpha_torus = alpha;
revs_per_second_torus = revs_per_sec;
r0_torus = alpha;           % assumes outer radius == 1

period_torus = period;

xc0_torus = 1;    % center of sphere used to initialize problem
yc0_torus = 0;
zc0_torus = alpha_torus;

if (T > 0)
    % Trace each point back to the origin
    reverse_torus = true;
    
    [xcm,ycm] = meshgrid(xcenter,ycenter);    
    x0 = xcm(:);
    y0 = ycm(:);   
    
    f0 = [x0; y0];   % Computational coordinates
    opt.RelTol = 1e-12;
    [~,pout] = ode45(g_psi_rhs,[0,T],f0,opt);
    
    N = length(x0);
    xc = pout(end,1:N);
    yc = pout(end,(N+1):end);
    
    mx = length(xcenter);
    my = length(ycenter);
    xcm = reshape(xc,mx,my);
    ycm = reshape(yc,mx,my);
    
    [xp,yp,zp] = g_mapc2m(xcm,ycm,alpha_torus);
    
    q = solution_init(xp,yp,zp);   % Computational boundaries        
else
    % Return initial conditions
    [xcm,ycm] = meshgrid(xcenter,ycenter);
    [xp,yp,zp] = g_mapc2m(xcm,ycm,alpha_torus);
    q = solution_init(xp,yp,zp);
end

if (nargout > 0)
    % We are plotting the solution on top of a ForestClaw solution
    qtrue = q;
else
    % This is being called as a standalone function.
    plot_torus_soln(xedge,yedge,q);
    hold on    
end

end

% -----------------------------------------------
% Plot the torus solution
% -----------------------------------------------
function plot_torus_soln(xedge,yedge,q)

global alpha_torus g_mapc2m

[xem,yem] = meshgrid(xedge,yedge);
[xep,yep,zep] = g_mapc2m(xem,yem,alpha_torus);

p = patch(surf2patch(xep,yep,zep));

set(p,'cdata',q);
set(p,'CDataMapping','scaled');  % Scale into current color map.
set(p,'FaceColor','flat');       % Single color per cell
set(p,'edgecolor','none');

daspect([1,1,1])
% set(p,'tag','torus');

end

% -----------------------------------------------
% Solution at the initial time
% -----------------------------------------------
function q = solution_init(x,y,z)

global alpha_torus xc0_torus yc0_torus zc0_torus;

x0 = xc0_torus;
y0 = yc0_torus;
z0 = zc0_torus;

r0 = alpha_torus;
Hsmooth = @(r) (tanh(r/0.02) + 1)/2;

r = sqrt((x - x0).^2 + (y - y0).^2 + (z - z0).^2);
q = Hsmooth(r + r0) - Hsmooth(r - r0);

end


%{
% -----------------------------------------------
% Get an initial set of x values. 
% -----------------------------------------------
function [xdata,ydata] = torus_init(xcenter,ycenter)

[xcm,ycm] = meshgrid(xcenter,ycenter);

xdata = xcm(:);
ydata = ycm(:);

end
%}

% -----------------------------------------------------------------
% Compute RHS : 
%
%      X' = u(X,Y)
%      Y' = v(X,Y)
% 
% Compute (u,v) from the stream function
% 
%      psi = (pi2*revs_per_s)*alpha*(pi2*eta + alpha*sin(pi2*eta))
% 
% -----------------------------------------------------------------
function fp = psi_rhs_torus(t,sigma)

global alpha_torus revs_per_second_torus reverse_torus;

revs_per_s = revs_per_second_torus;
alpha = alpha_torus;

N = length(sigma)/2; 

xc1 = sigma(1:N)';
yc1 = sigma(N+1:end)';

pi2 = 2*pi;

% Gradient psi in computational coordinates
psi_xi = 0;
psi_eta = (pi2)^2*revs_per_s*alpha*(1 + alpha*cos(pi2*yc1));

R = 1 + alpha*cos(pi2*yc1);
Reta = -pi2*alpha*sin(pi2*yc1);

% # T_xi
tau1 = zeros(3,N);
tau1(1,:) = -pi2*R.*sin(pi2*xc1);
tau1(2,:) =  pi2*R.*cos(pi2*xc1);
tau1(3,:) = 0;

% # T_eta
tau2 = zeros(3,N);
tau2(1,:) = Reta.*cos(pi2*xc1);
tau2(2,:) = Reta.*sin(pi2*xc1);
tau2(3,:) = pi2*alpha*cos(pi2*yc1);

% Compute u dot grad q in computational coordinates
c = cross(tau1,tau2);
w = 1./sqrt(dot(c,c));

u = w.*psi_eta;
v = -w.*psi_xi;

fp = sigma;
fp(1:N) = u(:);
fp(N+1:end) = v(:);

if reverse_torus
    fp = -fp;   
end

end

% -----------------------------------------------------------------
% Compute RHS : 
%
%      X' = u(X,Y)
%      Y' = v(X,Y)
% 
% Compute (u,v) from the stream function
% 
%    psi = pi2*revs_per_s*alpha*(pi2*(xc1+yc1) + alpha*sin(pi2*(xc1+yc1)))
% 
% -----------------------------------------------------------------
function fp = psi_rhs_twisted_torus(t,sigma)

global alpha_torus revs_per_second_torus reverse_torus;

revs_per_s = revs_per_second_torus;
alpha = alpha_torus;

N = length(sigma)/2; 

xc1 = sigma(1:N)';
yc1 = sigma(N+1:end)';

pi2 = 2*pi;


psi_xi = pi2^2*revs_per_s*alpha*(1 + alpha*cos(pi2*(xc1+yc1)));
psi_eta = pi2^2*revs_per_s*alpha*(1 + alpha*cos(pi2*(xc1+yc1)));

%         # Coordinate normals
R    = 1 +  alpha*cos(pi2*(xc1 + yc1));
Rxi  = -pi2*alpha*sin(pi2*(xc1 + yc1));
Reta = -pi2*alpha*sin(pi2*(xc1 + yc1));

%         # T_xi
tau1(1,:) = Rxi.*cos(pi2*xc1) - pi2*R.*sin(pi2*xc1);
tau1(2,:) = Rxi.*sin(pi2*xc1) + pi2*R.*cos(pi2*xc1);
tau1(3,:) = pi2*alpha*cos(pi2*(xc1 + yc1));

%         # T_eta
tau2(1,:) = Reta.*cos(pi2*xc1);
tau2(2,:) = Reta.*sin(pi2*xc1);
tau2(3,:) = pi2*alpha*cos(pi2*(xc1+yc1));

% Compute u dot grad q in computational coordinates
c = cross(tau1,tau2);
w = 1./sqrt(dot(c,c));

u = w.*psi_eta;
v = -w.*psi_xi;

fp = sigma;
fp(1:N) = u(:);
fp(N+1:end) = v(:);

if reverse_torus
    fp = -fp;
end

end



