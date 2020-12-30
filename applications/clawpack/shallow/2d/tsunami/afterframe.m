g = 1;
h0 = 1;
b = 12.2;
L = 4096;

plot_swe = true;

% set axis limits
if t <= 50
    w = 50;
    xticks = -w:10:w;
elseif (t < 75)
    w = 100;
    xticks = -w:20:w;
else
    w = 200;
    xticks = -10:50:w;
end

x1 = xticks(1);
x2 = xticks(end);
if (mq == 1)
    % axis([-w,w, -0.02, 0.1]);    
    axis([x1,x2, -0.02, 0.15]);
elseif mq == 2
    axis([x1,x2,-3.5,1]);
end

xl = [x1,x2];
yl = [-0.02,0.15];

% Set up moving reference frame
c = sqrt(g/h0);
refframe = c*t;

% Add plots
clear hvec lstr;
if (PlotType == 1)
    % axis([refframe-w/2, refframe+w/2,0,16]);
    set(gca,'xlim',[x1,x2]);
    ca = yl;
    caxis(ca);
    cv = linspace(0,ca(end),15);
    drawcontourlines(cv);
    colorbar;
    view(2)
elseif (PlotType == 4)
    hvec(1) = getlegendinfo();  % Get info from basic plot

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
                
    if (plot_swe)
        [amrdata_swe,t_swe] = readamrdata(1,Frame,'sgn1d/_output', ...
            'ascii');
        tol = 1e-6;
        if (abs(t_swe - t) > tol)
            fprintf('--------> WARNING : Times do not match; t_swe= %g; t = %g\n', t_swe, t);
        end
        hold on;
        if (UserVariable == 1)
            [x_swe,~,h_swe] = plotframe1ez(amrdata_swe,mq,'r-',UserVariableFile);
        else
            [x_swe,~,h_swe] = plotframe1ez(amrdata_swe,mq,'r-');
        end
        set(hvec,'Tag','SWE');
        hvec(2) = h_swe;
    end
    
    % hvec = get(gca,'children');
    
    for i = 1:length(hvec)
        h = hvec(i);
        if (~ishandle(h))
            continue;
        end
        tag = get(h,'Tag');
        if (strcmp(tag,'SWE') == 0)
            % set(h,'color','k')
        end
        set(h,'linewidth',2);
        xdata = get(h,'xdata');
        % set(h,'xdata',xdata/h0 - refframe);
        xlim([-10,100]);
    end
    
       
    % axis([-w,w,ylim]);
    set(gca,'xtick',xticks);
    ylim(yl);
        
    
    
    % mask = [1,plot_swe, plot_bouss, plot_basilisk] == 1;
    mask = hvec ~= 0;
    
    level = round(log2(mx));
    lstr{1} = sprintf('2d results (level %d)',level);

    if (plot_swe)
        mx_swe = amrdata_swe.mx;
        level_swe = round(log2(mx_swe));
        lstr{2} = sprintf('SWE (level %d)',level_swe);
    end
    set(gcf,'position',[360   509   835   189]);        
end

    
fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);


tstr = sprintf('t = %12.4f',t);
title(tstr);

set(gca,'fontsize',16);
set(gca,'box','on');

legend(hvec,lstr,'fontsize',16);

prt = false;
if (prt)
    fname = sprintf('plot_%04d.png',t);
    print('-dpng','-r512',fname);
end

shg

hold off;
