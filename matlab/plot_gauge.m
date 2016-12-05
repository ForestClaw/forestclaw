function hout = plot_gauge(gnum,mq)

if (nargin < 2)
    mq = 4;
end

l = {'Height (m)','u-momentum','v-momentum','Sea surface height (eta)'};

gdata = read_gauge_data(gnum);

t = gdata(:,2);
eta = gdata(:,mq+2);

h = plot(t,eta,'k');
xlabel('Time (s)','fontsize',16);
ylabel(l{mq},'fontsize',16);
title(sprintf('Gauge %d',gnum),'fontsize',18);
set(gca,'fontsize',16);
shg

if (nargout > 0)
    hout = h;
end

end
