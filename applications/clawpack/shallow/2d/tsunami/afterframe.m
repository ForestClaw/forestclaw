g = 1;
h0 = 1;
b = 12.2;

L = 4000;

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
        [x_swe,~,h_swe] = plotframe1ez(amrdata_swe,mq,'b-',UserVariableFile);
    else
        [x_swe,~,h_swe] = plotframe1ez(amrdata_swe,mq,'b-');
    end
    hvec(2) = h_swe;
end


if t <= 50
    w = 50;
    xticks = -w:10:w;
else
    w = 100;
    xticks = -w:20:w;
end


% mask = [1,plot_swe, plot_bouss, plot_basilisk] == 1;
mask = hvec ~= 0;

c = sqrt(g/h0);
refframe = c*t;

if (mq == 1)
    % axis([-w,w, -0.02, 0.1]);
    axis([-w,w, -0.02, 0.15]);
elseif mq == 2
    axis([-w,w,-3.5,1]);
end

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

set(gcf,'position',[360   509   835   189]);

set(gca,'fontsize',15);
tstr = sprintf('t = %12.4f',t);
title(tstr);
set(gca,'box','on');


% level = round(log(mx)/log(2));
% lstr{1} = sprintf('2d results (level %d)',level);
% 
% if (plot_swe)
%     lstr{2} = sprintf('SWE (level %d)',level_swe);
% end

prt = false;
if (prt)
    fname = sprintf('plot_%04d.png',t);
    print('-dpng','-r512',fname);
end

shg

hold off;
