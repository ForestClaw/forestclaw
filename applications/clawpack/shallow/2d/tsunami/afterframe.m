g = 1;
h0 = 1;
b = 12.2;
L = 4096;

plot_sgn1d = true;
plot_surf = true;

% set axis limits
if t <= 50
    w = 50;
    xticks = -w:10:w;
else
    w = 60;
    xticks = -w:20:w;
end

xl = [-w,w];
yl = [-0.02,0.1];

% Set up moving reference frame
c = sqrt(g/h0);
refframe = c*t;

shift_x = @(x) (x-b)/h0 - t*sqrt(g/h0);

% Add plots
if (PlotType == 1)
    % axis([refframe-w/2, refframe+w/2,0,16]);
    ca = yl;
    caxis(ca);
    cv = linspace(0,ca(end),15);
    cv(1) = [];
    drawcontourlines(cv);
    colorbar;
    set(gca,'xlim',[0,128]);
    set(gca,'ylim',[0,128]);
    
    showpatchborders;
    setpatchborderprops('linewidth',1);
    
    view(2)
    daspect([1,1,1]);
    
elseif (PlotType == 4)
    if (mq == 1)
        axis([xl,yl]);
    elseif mq == 2
        axis([-w,w,-3.5,1]);
    end

    [hvec,lstr] = getlegendinfo(0);  % Get info from basic plot
    set(hvec,'Tag','2d');

    ms = 12;  % Markersize
    userdata = get(gcf,'userdata');
    udlines = userdata.lines;
    maxlines = length(udlines);
    for i = 1:length(amrdata)
        level = amrdata(i).level;
        for j = 1:maxlines
            if (j-1 == level)
                set(udlines{j},'markersize',ms);  % Array
            end
        end
    end
                
    if (plot_sgn1d)
        [amrdata_sgn1d,t_sgn1d] = readamrdata(1,Frame,'sgn1d/_output', ...
            'ascii');
        tol = 1e-6;
        if (abs(t_sgn1d - t) > tol)
            fprintf('--------> WARNING : Times do not match; t_sgn1d = %g; t = %g\n', ...
                t_sgn1d, t);
        end
        hold on;
        if (UserVariable == 1)
            [q_sgn1d,x_sgn1d,h_sgn1d] = plotframe1ez(amrdata_sgn1d,mq,'r-',UserVariableFile);
        else
            [q_sgn1d,x_sgn1d,h_sgn1d] = plotframe1ez(amrdata_sgn1d,mq,'r-');
        end
        hvec(end+1) = h_sgn1d;
        set(hvec(end),'Tag','SGN1d');
        lstr{end+1} = sprintf('1d reference solution');
    end
    
    % hvec = get(gca,'children');
    
    set(gca,'xlim',[-20,128]);
    for i = 1:length(hvec)
        h = hvec(i);
        if (~ishandle(h))
            continue;
        end
        tag = get(h,'Tag');
        if (strcmp(tag,'SGN1d') ~= 0)
            % hvec contains all 1d data
            set(h,'linewidth',1,'color','k');
            xdata = get(h,'xdata');
            set(h,'xdata',shift_x(xdata));
        end 
    end
    % Set 2d line properties
    userdata = get(gcf,'userdata');
    udlines = userdata.lines;
    maxlevel = length(udlines);
    for level = 1:maxlevel
        level_data = udlines{level};
        if (~isempty(level_data))
            for k = 1:length(level_data)
                xdata = get(level_data(k),'xdata');
                % set(level_data(k),'color','g','marker','o','markersize',8);
                set(level_data(k),'xdata',shift_x(xdata));
            end
        end
    end
    
       
    % mask = [1,plot_sgn1d, plot_bouss, plot_basilisk] == 1;
    mask = hvec ~= 0;
    
    % lstr{1} = sprintf('2d results');

    if (plot_sgn1d)
    end
    
    axis([xl,yl]);
    set(gca,'xtick',xticks);
        
    set(gcf,'position',[2110,   509,   835,   189]);
    legend(hvec,lstr,'fontsize',12,'location','northwest');
        
    hold off;
    
end

    
fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

tstr = sprintf('t = %12.4f',t);
title(tstr);

set(gca,'fontsize',16);
set(gca,'box','on');

if plot_surf
    figure(2);
    for i = 1:length(amrdata)
        data = amrdata(i);
        xlow = data.xlow;
        ylow = data.ylow;
        dxp = data.dx;
        dyp = data.dy;
        xhi = xlow + dxp*mx;
        yhi = ylow + dyp*my;
        xe = linspace(xlow,xhi,mx+1);
        xc = xe(1:end-1) + dxp/2;
        ye = linspace(ylow,yhi,my+1);
        yc = ye(1:end-1) + dyp/2;
        
        [xm,ym] = meshgrid(xc,yc);
        qpatch = surface_eta(data.data(1,:)',xm,ym);
        qm = reshape(qpatch,size(xm));
        h = surf(shift_x(xm),ym,qm');
        set(h, 'edgecolor','none');
        hold on;        
    end
    
    ym = 64;
    plot3(shift_x(x_sgn1d),0*x_sgn1d + 64,q_sgn1d+1e-4,'r-','linewidth',2);
    set(gca,'xlim',xl);
    set(gca,'zlim',yl);
    set(gca,'xtick',xticks);
    daspect([10,50,0.05])
    view([-10,20]);
    hold off;
    figure(1);
end
    
prt = false;
if (prt)
    fname = sprintf('plot_%04d.png',t);
    print('-dpng','-r512',fname);
end

shg

hold off;
