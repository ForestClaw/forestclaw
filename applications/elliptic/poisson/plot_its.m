function plot_its(with_pc)


if (nargin < 1)
    close all;
    with_pc = false;
else
    hold on;
end

data = load('its.txt');

if (with_pc)
    c = 'r.-';
else
    c = 'b.-';
end

p = semilogy(data(:,1),data(:,2),c,'markersize',10);


title('Iterations','fontsize',18);
xlabel('N','fontsize',16);
ylabel('Residual','fontsize',16);
set(gca,'fontsize',16);



end