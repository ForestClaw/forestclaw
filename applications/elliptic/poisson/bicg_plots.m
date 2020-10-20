function bicg_plots()

% ------------------------ no preconditioning -----------------------------
figure(1);
clf
data_f = load('its_bicg_fortran.dat');

p(1) = semilogy(data_f(:,1),data_f(:,3),'r.-','markersize',10);
hold on;

data_te = load('its_bicg_te.dat');
p(2) = semilogy(data_te(:,1),data_te(:,2),'b.-','markersize',10);

data_te4x4 = load('its_bicg_te_4x4.dat');
p(3) = semilogy(data_te4x4(:,1),data_te4x4(:,2),'g.-','markersize',10);

semilogy([0 1050],[1,1]*1e-14,'k--');

lh = legend(p,{'Fortran (13s)','ThunderEgg (308s)','ThunderEgg 4x4 (210s)'},...
    'location','northeast');

xlim([0,1050]);

title('BICG - No preconditioner','fontsize',18);
xlabel('Iterations','fontsize',16);
ylabel('Residual','fontsize',16);
set(gca,'fontsize',16);
shg

% --------------------- with preconditioning -----------------------------
figure(2);
clf
data_f = load('its_bicg_pc_fortran.dat');

p(1) = semilogy(data_f(:,1),data_f(:,3),'r.-','markersize',10);
hold on;

% data_te = load('its_bicg_pc_te.dat');
% p(2) = semilogy(data_te(:,1),data_te(:,2),'b.-','markersize',10);

data_te4x4 = load('its_bicg_pc_te_4x4.dat');
p(2) = semilogy(data_te4x4(:,1),data_te4x4(:,2),'g.-','markersize',10);

semilogy([0 100],[1,1]*1e-14,'k--');

lh = legend(p,{'Fortran (4s)','ThunderEgg 4x4 (3600s)'},...
    'location','southeast');

xlim([0,100]);

title('BICG - With preconditioning','fontsize',18);
xlabel('Iterations','fontsize',16);
ylabel('Residual','fontsize',16);
set(gca,'fontsize',16);
shg


end