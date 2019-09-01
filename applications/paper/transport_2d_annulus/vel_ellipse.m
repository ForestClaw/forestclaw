function f = vel_ellipse(t,u)

global theta beta

pi2 = 2*pi;

[example,A,~,~,beta,theta,~,freq,cart_speed] = read_vars();


tfinal = 0.25;

N = length(u)/2;

f = zeros(size(u));

if example == 0
    % Rigid body rotation
    u1 = -1;  % Revs per second
    u2 = 0;
    X = u(1:N);
    Y = u(N+1:end);
    t1 = theta(1);
    t2 = theta(2);
    t0 = atan2(Y,X);
    m = t0 < 0;
    t0(m) = t0(m) + 2*pi;
    t0 = t0/(2*pi);
    xm = (t0-t1)/(t2-t1);
    rm = sqrt(X.^2 + Y.^2);
    ym = (rm-beta)/(1-beta);
    
    [t11,t21,t12,t22] = arrayfun(@covariant_basis,xm,ym);
    vc1 = u1.*t11 + u2.*t12;
    vc2 = u1.*t21 + u2.*t22;
else
    if example == 1
        vc1 = cart_speed;
        vc2 = 0;                
    elseif example == 2
        vc1 = cart_speed;
        vc2 = pi2*A*cos(freq*pi2*t/tfinal)/tfinal;
    elseif example == 3
        vc1 = cart_speed*pi*sin(pi*t/tfinal)/2.d0;
        vc2 = 0;
    elseif example == 4
        % Vertical motion
        vc1 = pi2*A*cos(freq*pi2*t/tfinal)/tfinal;
%        vc2 = cart_speed*pi*sin(pi*t/tfinal)/2.d0;
        vc2 = cart_speed;
    elseif example == 5
        X = u(1:N);
        Y = u(N+1:end);
        
        t1 = theta(1);
        t2 = theta(2);
        t0 = atan2(Y,X);
        m = t0 < 0;
        t0(m) = t0(m) + 2*pi;
        t0 = t0/(2*pi);
        xm = (t0-t1)/(t2-t1);
        rm = sqrt(X.^2 + Y.^2);
        ym = (rm-beta)/(1-beta);
                
        th = theta(1) + (theta(2) - theta(1))*xm;
        r = beta + (1-beta)*ym;
        w = 0.5;
        
        vc1 = -(1-w)*pi2*r.*sin(pi2*th);
        vc2 =  (1+w)*pi2*r.*cos(pi2*th);
        nc = sqrt(vc1.^2 + vc2.^2);
        
        vc1 = -vc1./nc;
        vc2 = -vc2./nc;   
    end
end

f(1:N) = vc1;
f(N+1:end) = vc2;
    
end

function [t11, t21, t12, t22] = covariant_basis(x,y)

global theta beta

pi2 = 2*pi;

tv = pi2*(theta(1) +  (theta(2)-theta(1))*x);
thetax = pi2*(theta(2)-theta(1));

r = beta + (1-beta)*y;
ry = 1-beta;

t11 = -thetax*r*sin(tv);
t21 = thetax*r*cos(tv);
    
t12 = ry*cos(tv);
t22 = ry*sin(tv);

end
    
