function plot_mixing
close all


fname = cell(6,1);

fname{1} = 'mix-CLAW-dx1.5-sp-CFL0.95.dat';
fname{2} = 'mix-CLAW-dx0.75-sp-CFL0.95.dat';
fname{3} = 'mix-CLAW-dx0.28125-sp-CFL0.95.dat';
fname{4} = 'mix-CLAW-dx1.5-un-CFL0.95.dat';
fname{5} = 'mix-CLAW-dx0.75-un-CFL0.95.dat';
fname{6} = 'mix-CLAW-dx0.1875-un-CFL0.95.dat';

fname{1} = 'mix.out';

iframe = load('diag.dat');

a = -0.8;
b = 0.9;

for i = 1:1,
  figure(i);
  clf;
  d = load(fname{i});
  xk = d(:,1);
  yk = d(:,2);
  area = d(:,3);
  [Lr,Lu,Lo,mr,mu,mo,sumA] = compute_diag(xk,yk,a,b,area);

  if (sumA < 1e-13)
    Lr = sumA;
    Lu = 0;
    Lo = 0;
  end

  fprintf('\n\n\n');
  fprintf('--------------------------------------------------------\n');
  fprintf('File : %s\n',fname{i});

  plot(xk(mr),yk(mr),'r.');
  hold on;
  plot(xk(mu),yk(mu),'g.');
  plot(xk(mo),yk(mo),'b.');

  fprintf('                                 Diagnostic     Fraction\n');
  fprintf('--------------------------------------------------------\n');
  fprintf('%30s %12.4e %12.4f\n','Real mixing (r)',Lr,Lr/sumA);
  fprintf('%30s %12.4e %12.4f\n','Range preserving unmixing (g)',Lu,Lu/sumA);
  fprintf('%30s %12.4e %12.4f\n','Over- and under-shooting (b)',Lo,Lo/sumA);
  fprintf('--------------------------------------------------------\n');
  fprintf('%30s %12.4e %12.4f\n','Total',sumA,1.0);

  % cosine bells are in range [Xmin, Xmax]
  Xmin = 0.1;
  Xmax = 1;

  % Correlated cosine bells are in range [ximin, ximax]
  ximin = a*Xmax^2 + b;
  ximax = a*Xmin^2 + b;

  r = linspace(0.1,1,201)';
  rl = (ximax - ximin)/(Xmin - Xmax)*(r - Xmin) + ximax;
  plot(r,a*r.^2 + b,'k-','linewidth',3)
  hold on;
  plot(r,rl,'k-','linewidth',3);
  set(gca,'xtick',0:0.1:1);
  set(gca,'ytick',0:0.1:1);

  xgmin = -0.15;
  xgmax = 1.1;
  ygmin = 0;
  ygmax = 1.1;
  axis([xgmin xgmax ygmin ygmax]);

  xlim = get(gca,'xlim');
  ylim = get(gca,'ylim');

  plot([Xmin Xmin],ylim,'k--');
  plot([Xmax Xmax],ylim,'k--');
  plot(xlim,[ximin ximin],'k--');
  plot(xlim,[ximax ximax],'k--');

  daspect([1 1 1]);

  xlabel('\phi','fontsize',14);
  ylabel('\Psi(\phi)','fontsize',14);
  title(sprintf('iframe = %d',iframe),'fontsize',18);

  prt = false;
  if (prt == true)
    filename = [fname{i} '.png'];
    prtstr = ['print -dpng ' filename];
    fprintf('%s\n',prtstr);
    print(filename,'png');
  end
  fprintf('--------------------------------------------------------\n');

end

end

% ------------------------------------------------------------------
function [Lr,Lu,Lo,mr,mu,mo,sumA] = compute_diag(xk,yk,a,b,area)

Xmin = 0.1;
Xmax = 1;
ximin = a*Xmax^2 + b;
ximax = a*Xmin^2 + b;

Rx = (Xmax - Xmin);
Ry = (ximax - ximin);

lk = ((ximax - ximin)/(Xmin - Xmax))*(xk - Xmin) + ximax;
curve_k = a*xk.^2 + b;

in_box = (xk >= Xmin) & (xk <= Xmax) & (yk >= ximin) & (yk <= ximax);
mr = in_box & ((yk >= lk) & (yk <= curve_k));
mu = in_box & ((yk > curve_k) | (yk < lk));
mo = (~mr) & (~mu);
mo = ~in_box;


% Now compute minimum distance.
c1 = 1/60;
c2 = 65340;
c3 = 29648025;

ck = c1*(c2*xk + 12*sqrt(12*(125*yk - 52).^3 + c3*xk.^2)).^(1/3);
x_root = ck + (1./ck).*(13/75 - (5/12)*yk);

x_root(Xmin > x_root) = Xmin;
x_root(Xmax < x_root) = Xmax;
y_root = a*x_root.^2 + b;

dk = sqrt(((xk - x_root)/Rx).^2 + ((yk - y_root)/Ry).^2);

A = sum(area);
Lr = sum(dk(mr).*area(mr))/A;
Lu = sum(dk(mu).*area(mu))/A;
Lo = sum(dk(mo).*area(mo))/A;
sumA = sum(dk.*area)/A;

K = length(xk);
if ((sum(mu) + sum(mr) + sum(mo)) ~= K)
  fprintf('\n');
  fprintf('sum(mr)   = %d\n',sum(mr));
  fprintf('sum(mu)   = %d\n',sum(mu));
  fprintf('sum(mo)   = %d\n',sum(mo));
  fprintf('Total     = %d\n',sum(mr) + sum(mu) + sum(mo));
  fprintf('K         = %d\n',K);
  error('compute_diag : Not all points accounted for');
end

end
