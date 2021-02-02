function alpha_conv(e)

close all;
hold on;

R_cyl = 100;
H_cyl = 2*pi;  % so that beta^2/4 = 1

beta = H_cyl/R_cyl;
if (nargin == 0)
    c = beta^2/4;
    e = sqrt(eps(c)/10);
end

fprintf('%25s %35.16e\n','e',e);
fprintf('%25s %35.16e\n','c + e^2',c + e^2);
fprintf('%25s %35.16e\n','c',c);


% Check
if c + e^2 ~= c
    fprintf('----------> e is not small enough');
    fprintf('\n');
end

% Need e <= 1e-8;
xbar = R_cyl*(1 + e);

alpha = xbar/R_cyl;
beta = H_cyl/R_cyl;

R_sphere = -R_cyl/2*(1 - alpha + beta^2/(4*(1-alpha)));
L = -R_cyl/2*(1 + alpha + beta^2/(4*(1-alpha)));



fprintf('%30s %30.16f\n','',)



%{

% Test for theta = 0 : 
ebar = -(1 - alpha);
w = ebar-e;
c = beta^2/4;
fprintf('%25s %35.16f\n','c',c);

fprintf('Approximation to R_sphere\n');
fprintf('%25s %30.16f\n','0.5*(e + w + c/e*(1-w/e))',0.5*(e + w   + c/e*(1-w/e)));  % Rs
fprintf('%25s %30.16f\n','0.5*(e + c/e - c*w/e^2)',  0.5*(e + c/e - c*w/e^2));  % Rs
fprintf('%25s %30.16f\n','0.5*(e + c/e - c*w/e^2)',  0.5*(    c/e - c*w/e^2));  % Rs
fprintf('%25s %30.16f\n','0.5*c/e*(1 - w/e)',        0.5*c/e*(1 - w/e));  % Rs
Rapprox = 0.5*c/e*(1 - w/e);
fprintf('\n');


fprintf('Approximation to L\n');
fprintf('%35s %30.16f\n','-(2 + e - c/e*(1-w/e))',-0.5*(2 + e - c/e*(1-w/e)));
fprintf('%35s %30.16f\n','-0.5*c/e*(2/c*e - 1 + c*w/e^2))',-0.5*c/e*(2/c*e  - 1 + w/e));
Lapprox = -0.5*c/e*(2/c*e  - (1 - w/e));
fprintf('\n');

fprintf('%25s %30.16f\n','Rs (approx)',Rapprox);
fprintf('%25s %30.16f\n','-L (approx)',-Lapprox);
fprintf('%25s %30.16f\n','Rs-L (approx)',Rapprox-Lapprox);
fprintf('\n');


fprintf('Need (R_sphere) : \n');
fprintf('%25s %30.16f\n','(e + c/e - c*w/e^2)',(e + c/e - c*w/e^2));
fprintf('%25s %30.16f\n','(e + c/e - c*w/e^2)',(c/e - c*w/e^2));
fprintf('\n');

fprintf('Need (-L) : \n');
fprintf('%25s %30.16f\n','2 + e - c/e*(1-w/e)',2 + e - c/e*(1-w/e));
fprintf('%25s %30.16f\n','2 - c/e*(1-w/e)',2 - c/e*(1-w/e));
fprintf('\n');

fprintf('cw/e^2 = %30.16f\n',c*w/e^2);


fprintf('\n');
fprintf('%25s %30.16f\n','Rs',R_sphere);
fprintf('%25s %30.16f\n','-L',-L);
fprintf('%25s %30.16f\n','Rs-L',R_sphere-L);
%}






% Get range for phi : 
if (-L <= R_cyl)
    phi0 = asin(-H_cyl/(2*R_sphere));
    phi1 = -phi0;
else
    phi1 = acos((L+R_cyl)/R_sphere);
    phi0 = -phi1;
end
phi = linspace(phi0,phi1,201);


fprintf('\n');
fprintf('%25s %35.16e\n','phi0',phi0);
fprintf('%25s %35.16e\n','phi1',phi1);

%{
fprintf('\n');
H2 = H_cyl/2;
fprintf('%25s %30.16e\n','asin(H/(2*Rs))',asin(H2/R_sphere));
fprintf('%25s %30.16e\n','2*e/Rc*(1+w/e)',H2*2*e/(R_cyl*c)*(1 + w/e));
p0 = H2*2*e/(R_cyl*c)*(1 + w/e);
p1 = H2*2*e/(R_cyl*c)*(1);
fprintf('%25s %30.16e\n','cos(p0)',cos(p0));
fprintf('%25s %30.16e\n','cos(p1)',cos(p1));
%}


R = R_sphere*cos(phi) - L;

% fprintf('%20s  %24.16f\n','phi0',phi0);
% fprintf('%20s  %24.16f\n','phi1',phi1);
% fprintf('\n');
% fprintf('%20s  %24.17g\n','R_sphere',R_sphere);
% fprintf('%20s  %24.17g\n','-L (xcenter)',-L);
% fprintf('\n');

fprintf('\n');
fprintf('%25s %35.16f\n','min(R)',min(R));
fprintf('%25s %35.16f\n','max(R)',max(R));
fprintf('%25s %35.8e\n','rel. diff.',(max(R)-min(R))/max(R));
% fprintf('%20s %24.16f\n','(L+Rc)/Rs',(L+R_cyl)/R_sphere);
% b2 = beta^2;
% X = (b2 - 4*(1-alpha).^2)./(b2 + 4*(1-alpha).^2);
% fprintf('%20s %24.16f\n','X',X);
% 
fprintf('\n');


avec = linspace(-10,10,100);
R_sphere = -R_cyl/2*(1 - avec + beta^2./(4*(1-avec)));
L = -R_cyl/2*(1 + avec + beta^2./(4*(1-avec)));
plot(avec,(L + R_cyl)./R_sphere,'k');
hold on;
X = (1 - 4*((1-avec)/beta).^2)./(1 + 4*((1-avec)/beta).^2);
plot(avec,X,'b.','markersize',15);

xlabel('alpha','fontsize',16);
ylabel('cos(phi0)','fontsize',16);

% plot(phi,(L+R_cyl))

shg



end