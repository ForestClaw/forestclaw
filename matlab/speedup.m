% Compute speedup and efficiency
% 
% Model : S = T(1)/T(p) = p^alpha
% 
%   T(p) = T(1)*p^(-alpha)
%   
%       T = C*p^(-alpha)
%   log(T) = alpha*(-log(p)) + log(C)
%   
% We expect alpha \approx 1.
%

function alphaout = speedup(p,T)
figure(2);
clf;

T = T(:);
p = p(:);

a = polyfit(-log(p),log(T),1);

alpha = a(1);
C = exp(a(2));

loglog(p,T,'.','markersize',30);
hold on;
ph(1) = loglog(p,C*p.^(-alpha),'r');

ph(2) = loglog(p,T(1)./p,'k--');
loglog(p,T(1)./p,'k.','markersize',30);

title('Strong scaling','fontsize',18);
str1 = sprintf('Actual (\\alpha = %6.3f)\n',alpha);
str2 = sprintf('Optimal (\\alpha = 1)');
lh = legend(ph,{str1,str2},'location','southwest');
set(lh,'fontsize',16);


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

