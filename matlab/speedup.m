% SPEEDUP computes speedup and efficiency of a parallal process
%
% ALPHA=SPEEDUP(P,T) computes the best-fit speed-up and parallel efficiency
% of a process on with timings in vector T, on P processors.  The model
% used is : 
% 
%   Speedup :     S = T(1)/T(p) = p^alpha
%       -->    T(p) = T(1)*p^(-alpha)
%   
%   T = C*p^(-alpha)
%   log(T) = alpha*(-log(p)) + log(C)
%   
% We expect alpha \approx 1.
%

function alphaout = speedup(p,T)
T = T(:);
p = p(:);

a = polyfit(-log(p),log(T),1);

alpha = a(1);
C = exp(a(2));

figure(2);
clf;
loglog(p,T,'.','markersize',30);
hold on;
ph(1) = loglog(p,C*p.^(-alpha),'r');

ph(2) = loglog(p,T(1)./p,'k--');
loglog(p,T(1)./p,'k.','markersize',30);

title('Strong scaling (slope = \alpha)','fontsize',18);
str1 = sprintf('Actual (\\alpha = %6.3f)\n',alpha);
str2 = sprintf('Optimal (\\alpha = 1)');
lh = legend(ph,{str1,str2},'location','southwest');
xlabel('Proc count','fontsize',16);
ylabel('Speed-up','fontsize',16);
set(lh,'fontsize',16);
set(gca,'fontsize',16);

% efficiency vs. speedup.
clear ph;
figure(3)
clf
pv = logspace(0,5,1000);
T_llsq = C*pv.^(-alpha);
S_llsq = T(1)./T_llsq;
S_computed = T(1)./T;
S_true = pv;

eff_llsq = S_llsq./pv;
eff_computed = S_computed./p;
eff_true = S_true./pv;
semilogx(p,eff_computed,'.','markersize',30);
hold on;
ph(1) = semilogx(pv,eff_llsq,'b');
ph(2) = semilogx(pv,eff_true,'b--');
ylim([0 1.2]);


title('Efficiency (S/p)','fontsize',18);
str1 = sprintf('Actual');
str2 = sprintf('Optimal (eff. = 1)');
lh = legend(ph([1 2]),{str1,str2},'location','southwest');
set(lh,'fontsize',16);
xlabel('Proc count','fontsize',16);
ylabel('Efficiency','fontsize',16);
set(gca,'fontsize',16);


S = T(1)./T;
e = T(1)./T./p;

fprintf('\n');
fprintf('%2s %12s %16s %12s\n','p','T(p)','S(p)=T(1)/T(p)','e=S(p)/p');
fprintf('---------------------------------------------\n');
fprintf('%2d %12.2f %16.2f %12.3f\n',[p'; T'; S'; e']);
fprintf('\n');
fprintf('%s %12.4f\n','alpha',alpha);

if (nargout > 0)
    alphaout = alpha;
end

