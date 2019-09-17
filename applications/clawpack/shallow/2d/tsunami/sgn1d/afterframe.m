g = 1;
h0 = 1;
b = 12.2;

L = 4000;

hvec = zeros(1,4);
hvec(1) = getlegendinfo();  % Get info from basic plot

plot_swe = true;   % High res - 2^15

if (plot_swe)
    [amrdata_swe,t_swe] = readamrdata(1,Frame,'./_output_swe', ...
        'ascii');
    tol = 1e-6;
    if (abs(t_swe - t) > tol)
        error('Times do not match; ');
    end
    hold on;
    if (UserVariable == 1)
        [x_swe,~,h_swe] = plotframe1ez(amrdata_swe,mq,'b-',UserVariableFile);
    else
        [x_swe,~,h_swe] = plotframe1ez(amrdata_swe,mq,'b-');
    end
    hvec(2) = h_swe;
    level_swe = round(log(length(x_swe))/log(2));
end


c = sqrt(g/h0);
refframe = c*t;

if t <= 50
    w = 50;
    xticks = -w:10:w;
else
    w = 100;
    xticks = -w:20:w;
end


% mask = [1,plot_swe, plot_bouss, plot_basilisk] == 1;
mask = hvec ~= 0;

for i = 1:length(hvec)
    h = hvec(i);
    if (mask(i))
        set(h,'linewidth',1,'markersize',10);
        xdata = get(h,'xdata');
        set(h,'xdata',xdata - b - t);
    end
end

if (mq == 1)
    axis([-w,w, -0.02, 0.1]);
    % axis([-w,w, -0.2, 0.2]);
elseif mq == 2
    axis([-w,w,-3.5,1]);
end
tstr = sprintf('t = %.0f',t);

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

set(gcf,'position',[360   509   835   189]);

set(gca,'fontsize',15);
tstr = sprintf('t = %.0f',t);
title(tstr);
set(gca,'box','on');


level = round(log(mx)/log(2));
lstr{1} = sprintf('AMRClaw (level %d)',level);

if (plot_swe)
    lstr{2} = sprintf('SWE (level %d)',level_swe);
end



if length(hvec) == 4
    if (t < 1000)
        %legend(lstr{mask},'location','northeast');
    else
        %legend(lstr{mask},'location','northwest');
    end
end
% legend(hvec(mask),lstr{mask},'location','northwest');


prt = true;
if (prt)
    fname = sprintf('plot_%04d.png',t);
    print('-dpng','-r512',fname);
end

shg

hold off;
