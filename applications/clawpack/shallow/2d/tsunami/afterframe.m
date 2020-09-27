g = 1;
h0 = 1;
b = 12.2;
L = 4000;

c = sqrt(g/h0);
refframe = c*t + b;

if t <= 50
    w = 50;
    xticks = -w:10:w;
    ylim = [-0.02, 0.1];
else
    w = 100;
    xticks = -w:20:w;
    ylim = [-0.02, 0.1];
end


if (PlotType == 1)
    axis([refframe-w/2, refframe+w/2,0,16]);
    ca = [-0.02,0.05];
    caxis(ca);
    cv = linspace(0,ca(end),15);
    drawcontourlines(cv);
    colorbar;
    view(2)
elseif (PlotType == 4)
    hvec = zeros(1,2);
    hvec(1) = getlegendinfo();  % Get info from basic plot
    
    plot_swe = true;   % High res - 2^15
    
    if (plot_swe)
        [amrdata_swe,t_swe] = readamrdata(1,Frame,'sgn1d/_output', ...
            'ascii');
        tol = 1e-6;
        if (abs(t_swe - t) > tol)
            fprintf('Times do not match; t_swe= %g; t = %g\n', t_swe, t);
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
    
    hvec = get(gca,'children');
    
    for i = 1:length(hvec)
        h = hvec(i);
        tag = get(h,'Tag');
        if (strcmp(tag,'SWE') == 0)
            % set(h,'color','k')
        end
        set(h,'linewidth',2);
        xdata = get(h,'xdata');
        set(h,'xdata',xdata/h0 - refframe);
    end
    
       
    axis([-w,w,ylim]);
    set(gca,'xtick',xticks);
    
    
    % mask = [1,plot_swe, plot_bouss, plot_basilisk] == 1;
    mask = hvec ~= 0;
    
    level = round(log(mx)/log(2));
    lstr{1} = sprintf('2d results (level %d)',level);

    if (plot_swe)
        lstr{2} = sprintf('SWE (level %d)',level_swe);
    end
    set(gcf,'position',[360   509   835   189]);    
end
    
    
fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);


set(gca,'fontsize',15);
tstr = sprintf('t = %12.4f',t);
title(tstr);
set(gca,'box','on');


prt = false;
if (prt)
    fname = sprintf('plot_%04d.png',t);
    print('-dpng','-r512',fname);
end

shg

hold off;
