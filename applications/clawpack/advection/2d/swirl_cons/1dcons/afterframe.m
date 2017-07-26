hold on;
dir = './1dadv/';
dim = 1;
[amrdata1d,t1d] = readamrdata(dim,Frame,dir);
[q1d,x1d] = plotframe1ez(amrdata1d,mq,'b.-');

if (abs(t1d-t) > 1e-12)
    fprintf('t1d = %16.8e\n',t1d);
    fprintf('t   = %16.8e\n',t);
    error('Reference solution time incorrect');
end

% Compute errors
pp = csape(x1d,q1d,'periodic');
qxc = ppval(pp,xcenter);
qerr = abs(qxc(:)-q(:));
qnorm(1) = sum(qerr)*dx;
qnorm(2) = sqrt(sum(qerr.^2*dx));
qnorm(3) = max(qerr);
fprintf('%5d %12.4e %12.4e %12.4e %12.4e %12.4e\n',mx,qnorm(1),...
    qnorm(2),qnorm(3),qmin,qmax);

% Plot zeros of velocity function
vel_case = 3;
yl = [-5,25];
switch vel_case
    case 1
        % Case 1 : u = cos(2*pi*x)  + 2
        u = @(x) cos(2*pi*x) + a;
        z = acos(-a)/(2*pi);
        plot([z,z],yl,'k--');
        plot(1-[z,z],yl,'k--');

    case 2
        % Case 2 : u = cos(2*pi*x) 
        u = @(x) cos(2*pi*x) + a;
        z = acos(-a)/(2*pi);
        plot([z,z],yl,'k--');
        plot(1-[z,z],yl,'k--');
        
    case 3
        y1 = @(x) 0.5*tanh((x-0.65)/0.025) + 1;
        y2 = @(x) 0.5*tanh((0.35-x)/0.04) + 1;
        u = @(x) y1(x) + y2(x) - 1.5;
    otherwise
        %u = @(x) x*nan;
end

xc = linspace(0,1,200);
plot(xc,u(xc),'r--');   % Velocity
plot(xcenter,0*xcenter,'k');

fprintf('qmin = %12.4e\n',qmin);

qm = min([max([q1d(:); qmax]),15]);
ylim([-2 max([2.1,1.1*qm]);])

set(gca,'fontsize',16);
set(gca,'box','on');

shg

hold off;

clear afterframe
