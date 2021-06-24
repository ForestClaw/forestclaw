[example,A,rinit,beta,theta,freq,cart_speed] = read_vars();

setslicecolor('none');
hold on;


ravg = (1+beta)/2;

% Plot path
% Handled by afterframe

th = 2*pi*linspace(0.25-1/32,0.25+1/32,200);
xdisk = ravg*cos(th);
ydisk = ravg*sin(th);
hold on;
plot(xdisk,ydisk,'k','linewidth',2);
%plot(xdisk([1,end]),ydisk([1 end]),'k','linewidth',2);


% Plot starting disk
th = linspace(0,2*pi,200);
R = 0.05;
xc = xpath(1);
yc = ypath(1);
x0 = xc + R*cos(th);
y0 = yc + R*sin(th);
% plot(x0,y0,'k','linewidth',3);
z0 = x0*0 - 0.05;
fill3(x0,y0,z0,[1,1,1]*0.7);

% Plot ending disk
xc = xpath(end);
yc = ypath(end);
x0 = xc + R*cos(th);
y0 = yc + R*sin(th);
% plot(x0,y0,'k','linewidth',3);
z0 = x0*0 - 0.05;
fill3(x0,y0,z0,[1,1,1]*0.7);

% Plot path
plot_path = false;
if (plot_path && example ~= 5)    
    hold on;
    
    N = 128;
    th0 = linspace(0,2*pi,N+1);
    x0 = rinit*cos(th) + pstart(1);
    y0 = rinit*sin(th) + pstart(2);

    Y0 = pstart(:);
    [tout,yout] = ode45(@vel_ellipse,[0,tfinal],Y0);
    xpath = yout(:,1)';
    ypath = yout(:,2)';

    % plot patch
    hold on;
    plot([pstart(1),pend(1)],[pstart(2),pend(2)],'k.','markersize',20);
    plot(xpath,ypath,'k','linewidth',2);
    
    if (t > 0)
        Y0 = [x0(:); y0(:)];    
        [tout,yout] = ode45(@vel_ellipse,[0,t],Y0);
        xth = yout(end,1:N+1)';
        yth = yout(end,N+2:end)';
    else
        xth = x0;
        yth = y0;
    end
        
    plot(xth,yth,'k','linewidth',2);
    hold off
    
    if (example == 4)
        hold on;        
        t1 = theta(1);
        t2 = theta(2);
        t0 = atan2(y0,x0);
        m = t0 < 0;
        t0(m) = t0(m) + 2*pi;
        t0 = t0/(2*pi);
        xm = (t0-t1)/(t2-t1);
        rm = sqrt(x0.^2 + y0.^2);
        ym = (rm-beta)/(1-beta);
        % hq = ellipse_vel(4,1.5,xm,ym);
        % set(hq,'linewidth',2,'color','k');
        hold off
    end
end


showgridlines

axis off
delete(get(gca,'title'));


dir = '/Users/calhoun/projects/ForestClaw/papers/sisc_2014/figs_annulus/tmp/';
% Annulus_1 : constant theta
%print('-dpdf',[dir,'annulus_1.pdf'])

% Annulus_2 : constant r
print('-dpdf',[dir,'annulus_2.pdf'])


plot_zoom = false;
if (plot_zoom)
    axis([-0.036,   0.036, 0.6985, 0.7005]);
    set(gca,'dataaspectratiomode','auto')
    setslicecolor([1,1,1]*0.8);
    print('-dpdf',[dir,'gap.pdf'])
end

shg

    

