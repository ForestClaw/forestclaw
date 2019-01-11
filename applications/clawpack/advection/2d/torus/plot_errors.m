function plot_errors(e)
close all

Nvec = e(:,1);
n = length(Nvec);

c = {'r','b','g'};
idx = 1:n;
lh = zeros(3,1);
lstr = {'1-norm','2-norm','inf-norm'};
for i = 1:3
    p = polyfit(log(Nvec(idx)),log(e(idx,i+1)),1);
    lh(i) = loglog(Nvec,exp(polyval(p,log(Nvec))),c{i},'linewidth',2);  
    lstr{i} = sprintf('%10s (rate = %6.4f)',lstr{i},-p(1));
    hold on;
end
loglog(Nvec,e(:,2:4),'k.','markersize',30);
hold on;

i1 = idx(1);
i2 = idx(end);
yl = ylim;
plot([Nvec(i1) Nvec(i1)],yl,'k--');
plot([Nvec(i2) Nvec(i2)],yl,'k--');

legend(lh,lstr,'location','northeast');

xlabel('N','fontsize',16);
ylabel('Error','fontsize',16);
title('Error (Torus)','fontsize',18);

set(gca,'fontsize',16);

conv_rates = log(e(1:end-1,2:4)./e(2:end,2:4))/log(2);
fprintf('\n');
fprintf('          Convergence rates\n');
fprintf('%s\n',double('-')*ones(1,37));
cr = [Nvec(2:end) Nvec(1:end-1) conv_rates];
fprintf('%4d/%4d %8.4f %8.4f %8.4f\n',cr');


end