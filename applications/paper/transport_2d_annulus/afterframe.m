setviews;

% Geometry
d = load('_output/mapping.dat');
A = d(1);
rinit = d(2);
beta = d(3);
freq = d(6);

tfinal = 0.25;
vcart = 1.092505803290319;
dlen = vcart*tfinal;

ravg = (1 + beta)/2;
thc = pi/2*(1 + 1/8);
pstart = [ravg*cos(thc),ravg*sin(thc)];
pend = pstart + [dlen, 0];

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

% axis([-0.707106781186547   0.707106781186548   0.282842712474619,1]);

axis([-0.211745637207774, 0.211745637207774, 0.5, 0.85])


% showgridlines;

plot_path = true;
if (plot_path)    
    hold on;
    % Plot path
    tvec = linspace(0,tfinal,200);
    
    xpath = pstart(1) + tvec*vcart;
    ypath = pstart(2) + A*sin(2*pi*freq*tvec/tfinal);
    
    p_path = plot(xpath,ypath,'k','linewidth',2);
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
        
    xth = x0 + t*vcart;
    yth = y0 + A*sin(2*pi*freq*t/tfinal);
        
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
