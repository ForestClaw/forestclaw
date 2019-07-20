setviews;

% Geometry
d = load('_output/mapping.dat');
A = d(1);
rinit = d(2);
pstart = d(3:4)';
beta = d(5);
theta = d(6:7);
vcart = d(8:9);
freq = d(10);

tfinal = 0.25;
dlen = max(abs(vcart))*tfinal;

if (vcart(2) == 0)
    pend = pstart + [dlen, 0]*vcart(1)/abs(vcart(1));
else
    pend = pstart + [0,dlen]*vcart(2)/abs(vcart(2));
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

%axis([-0.707106781186547   0.707106781186548   0.282842712474619,1]);

rotate_position = 4*mean(theta) - 1;

if (rotate_position == 0)
    % Top
    axis([-0.211745637207774, 0.211745637207774, 0.5, 0.85])
elseif (rotate_position == 1)
    % left
    axis([-0.85,-0.5,-0.211745637207774, 0.211745637207774])
elseif (rotate_position == 2)
    % bottom
    axis([-0.211745637207774, 0.211745637207774, -0.85, -0.5])
else
    % right
    axis([0.5,0.85,-0.211745637207774, 0.211745637207774])
end
%axis image

% showgridlines;

plot_path = true;
if (plot_path)    
    hold on;
    % Plot path
    tvec = linspace(0,tfinal,200);
    
    %{
    if (vcart(2) == 0)
       xpath = pstart(1) + tvec*vcart(1);
       ypath = pstart(2) + A*sin(2*pi*freq*tvec/tfinal);         
    else
        xpath = pstart(1) + A*sin(2*pi*freq*tvec/tfinal);
        ypath = pstart(2) + tvec*vcart(2);
    end
    
    p_path = plot(xpath,ypath,'k','linewidth',2);
    %}
    hold on;

    % Plot start and endpoints
    plot(pstart(1),pstart(2),'k.','markersize',30);
    hold on;
    plot(pend(1),pend(2),'k.','markersize',30);
    plot([pstart(1),pend(1)],[pstart(2),pend(2)],'k','linewidth',2);
        
    N = 128;
    th = linspace(0,2*pi,N+1);
    x0 = rinit*cos(th) + pstart(1);
    y0 = rinit*sin(th) + pstart(2);

    if (vcart(2) == 0)
%         xth = x0 + t*vcart(1);
%         yth = y0 + A*sin(2*pi*freq*t/tfinal);
          xth = x0 + (dlen/2)*(1-cos(pi*t/tfinal));
          yth = y0;
    else
        xth = x0 + A*sin(2*pi*freq*t/tfinal);
        yth = y0 + t*vcart(2);
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
