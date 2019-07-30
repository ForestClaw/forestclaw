setviews;

% Geometry
[example,A,rinit,pstart, beta,theta,vcart, ...
        freq,cart_speed] = read_vars();

tfinal = 0.25;
vcart = cart_speed;
dlen = vcart*tfinal;

if (example < 4)
    pend = pstart(:)' + [dlen, 0];
else
    pend = pstart(:)' + [0, dlen];
end

fprintf('%10s %24.16e\n','qmin',qmin);
fprintf('%10s %24.16e\n','qmax',qmax);


s = 1e-2;    
alim = [-1,1];
axis([alim alim]);
daspect([1 1 1]);
view(vtop)

showpatchborders(1:10);
setpatchborderprops('linewidth',1)


caxis([-1,1]);

set(gca,'fontsize',16);

axis([-0.707106781186547   0.707106781186548   0.282842712474619,1]);

%axis([-0.211745637207774, 0.211745637207774, 0.5, 0.85])

plot_path = true;
if (plot_path)    
    hold on;
    
    N = 128;
    th0 = linspace(0,2*pi,N+1);
    x0 = rinit*cos(th0) + pstart(1);
    y0 = rinit*sin(th0) + pstart(2);
    
    Y0 = pstart(:);
    [tout,yout] = ode45(@vel_ellipse,[0,tfinal],Y0);
    xpath = yout(:,1)';
    ypath = yout(:,2)';

    % plot patch
    hold on;
    plot([pstart(1),pend(1)],[pstart(2),pend(2)],'k.','markersize',20);
    plot(xpath,ypath,'k','linewidth',2);

    
    if (t > 0)
        Y0 = [x0(:); y0(:);];    
        [tout,yout] = ode45(@vel_ellipse,[0,t],Y0);
        xth = yout(end,1:N+1)';
        yth = yout(end,N+2:end)';
    else
        xth = x0;
        yth = y0;
    end
        
    plot(xth,yth,'k','linewidth',2);
    hold off    
end


%
NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 8;
  axis([0 1 0 1]);
  filename = sprintf('annulus_%04d.png',Frame)
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
