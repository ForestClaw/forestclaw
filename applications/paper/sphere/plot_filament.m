clear all;
clf;

fname = cell(6,1);

fname{1} = 'lf-CLAW-dx1.5-sp-CFL0.95.dat';
lstr{1} = '\Delta\lambda = 1.50^o (sp)';

fname{2} = 'lf-CLAW-dx0.75-sp-CFL0.95.dat';
lstr{2} = '\Delta\lambda = 0.75^o (sp)';

fname{3} = 'lf-CLAW-dx0.28125-sp-CFL0.95.dat';
lstr{3} = '\Delta\lambda = 0.28^o (sp)';

fname{4} = 'lf-CLAW-dx1.5-un-CFL0.95.dat';
lstr{4} = '\Delta\lambda = 1.50^o (un)';

fname{5} = 'lf-CLAW-dx0.75-un-CFL0.95.dat';
lstr{5} = '\Delta\lambda = 0.75^o (un)';

fname{6} = 'lf-CLAW-dx0.1875-un-CFL0.95.dat';
lstr{6} = '\Delta\lambda = 0.19^o (un)';


s.sym = {'r.-','g.-','b.-','ro-','go-','bo-'};
s.size = [30, 30, 30, 10, 10, 10];

plot([0 1],[100 100],'k-','linewidth',3);
hold on;
plot([0 1],[0 0],'k-','linewidth',3);

idx = [1 4 2 5 3 6];
p = zeros(6,1);
for i = 1:6
  ii = idx(i);
  fn = fname{ii};
  d = load(fn);
  p(i) = ...
      plot(d(:,1),d(:,2),s.sym{ii},'markersize',s.size(ii),'linewidth',2);
  hold on;
end
set(gca,'xtick',[0:0.1:(1+1e-8)]);
axis([0.1 1 -5 130]);
grid on

lh = legend(p,lstr{idx},3);
set(lh,'fontsize',16);

title('Filament diagnostic (cosine hill)','fontsize',18);
xlabel('Filament threshold','fontsize',16);
set(gca,'fontsize',16);
