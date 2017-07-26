if PlotType == 1

    s = 0;
    axis([-s 1+s -s 1+s])
    daspect([1 1 1]);
    axis on;
    
    if (ShowUnderOverShoots == 1)
        qlo = 0;
        qhi = 1;
        under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
        over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
        fprintf('%-10s %12.4e\n','qmin',qmin);
        fprintf('%-10s %12.4e\n','qmax',qmax);
        colorbar_underover(under_label,over_label);
    else                 
        yrbcolormap;
        colorbar;
        caxis([0 1]);
    end
    showpatchborders;
    setpatchborderprops('linewidth',1);
    view(2);
elseif PlotType == 4
    hold on;
    dir = './1dcons/1dadv/';
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
    qerr = abs(qxc(:)-q(:,2));
    qnorm(1) = sum(qerr)*dx;
    qnorm(2) = sqrt(sum(qerr.^2*dx));
    qnorm(3) = max(qerr);
    fprintf('%5d %12.4e %12.4e %12.4e %12.4e %12.4e\n',mx,qnorm(1),...
        qnorm(2),qnorm(3),qmin,qmax);

    % Plot zeros of velocity function
    vel_case = 4;
    yl = [-5,25];
    switch vel_case
        case 1
            u = @(x) (x < 0.5)*(-0.5) + (x >= 0.5)*(0.5);
        case 2
            % Case 2 : u = cos(2*pi*x) + a
            a = 0;
            u = @(x) cos(2*pi*x) + a;
            z = acos(-a)/(2*pi);
            plot([z,z],yl,'k--');
            plot(1-[z,z],yl,'k--');
        case 3
            u = @(x) sin(2*pi*x).*sin(16*pi*x);
        case 4
            a = 0.01;
            u = @(x) 1 + 0.5*tanh((x-0.5)/a);
        case 5
            a = 0.01;
            u = @(x) tanh((x-0.5)/a);
        case 6
            a = 0.1;
            u = @(x) -tanh((x-0.5)/a);
        otherwise
            %u = @(x) x*nan;
    end
    
    xc = linspace(0,1,200);
    plot(xc,u(xc),'r--');   % Velocity
    plot(xcenter,0*xcenter,'k');    
    
    fprintf('qmin = %12.4e\n',qmin);
    
    qm = min([max([q1d(:); qmax]),15]);
    ylim([-2 max([2.1,1.1*qm]);])
       
    hold off;
     
end
    
set(gca,'fontsize',16);
set(gca,'box','on');

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'swirl0000','png');
  print('-dpng',filename);
end

shg

clear afterframe;
