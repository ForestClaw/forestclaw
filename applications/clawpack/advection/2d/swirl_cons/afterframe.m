if PlotType == 1

    s = 0;
    axis([-s 1+s -s 1+s])
    daspect([1 1 1]);
    axis on;
    
    hidepatchborders;
    
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
    showpatchborders(1:10);
    view(2);
elseif PlotType == 4
    hold on;
    dir = './1dadv/';
    dim = 1;
    [amrdata1d,t1d] = readamrdata(dim,Frame,dir);
    [q1d,x1d] = plotframe1ez(amrdata1d,mq,'b.-');    

    % Plot zeros of velocity function
    vel_case = 1;
    yl = [-5,25];
    switch vel_case
        case 1
            % Case 1 : u= cos(2*pi*x)
            u = @(x) cos(2*pi*x);
            plot([0.25,0.25],yl,'k--');
            plot([0.75,0.75],yl,'k--');
        case 2
            % Case 2 : u = cos(2*pi*x) + 0.5
            u = @(x) cos(2*pi*x) + 0.5;
            z = acos(-0.5)/(2*pi);
            plot([z,z],yl,'k--');
            plot(1-[z,z],yl,'k--');
        case 3
            % Case 3 : u = cos(2*pi*x) + 1
            u = @(x) cos(2*pi*x) + 1;
            plot([0.5,0.5],yl,'k--');
        case 4
            % Case 4 : u = cos(2*pi*x) + 2
            u = @(x) cos(2*pi*x) + 2;
            % no zeros
        case 5
            u = @(x) sin(2*pi*x).*sin(16*pi*x);
        otherwise
            u = @(x) (x < 0.5)*(-0.5) + (x >= 0.5)*(0.5);
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

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'swirl0000','png');
  print('-dpng',filename);
end

shg

clear afterframe;
