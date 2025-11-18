% Variable Coefficient Poisson Problem

function vc_poisson()

global b1 b2 u0_shape;

close all;

% Exact solution
ue = @(x,y) (1-x).*x.*(1-y).*y.*exp(x.*y);

% Beta
beta_coeffs = @(x,y) 1 + x.*y;

% Domain [a,b]x[a,b]
a = 0;
b = 1;

% Method choice 
use_cg = true;    % use_cg; use_bicg
use_bicg = ~use_cg;

np = 1;    % Norm to plot (1,2,3)

% Numerical parameters
N0 = 32;
nvec = N0*2.^(0:4);

tol = 1e-14;
maxit = 1e6;

area = (b-a)^2;

% Start convergence study

len = 50;
dashed_line = sprintf('%s',double('-')*ones(1,len));
fprintf('%s\n',dashed_line);
if (use_cg)
    fprintf('Conjugate Gradient\n');
else
    fprintf('BicgStab\n');
end
fprintf('%5s %10s %16s %16s\n','N','Iters.','Error','Time');
fprintf('%s\n',dashed_line);

for i = 1:length(nvec)
    N = nvec(i);   
    h = (b-a)/N;
    x = linspace(a,b,N+1);
    y = x;
    [xm,ym] = ndgrid(x,y);

    u_exact = ue(xm,ym);
    u_bc = u_exact;
    u_bc(2:N,2:N) = 0;

    xc = x(1:end-1) + h/2;
    [xcm,ycm] = ndgrid(xc,y);
    b1 = beta_coeffs(xcm,ycm);
    
    yc = y(1:end-1) + h/2;
    [xcm,ycm] = ndgrid(x,yc);
    b2 = beta_coeffs(xcm,ycm);

    f = frhs(xm,ym);
    
    %{
    figure(5);
    surf(xm,ym,bta);
    title('Variable coefficient beta');
    %}

    % For reshaping arrays
    u0_shape =  zeros(N+1,N+1);

    if use_cg 
        f = f - matvec(u_bc);
        f = reshape(f,size(u0_shape));
    elseif use_bicg
        f = f(:) - matvec(u_bc(:));
        f = reshape(f,size(u0_shape));
    end
    f([1 end],:) = 0;
    f(:,[1 end]) = 0;
    
    if (use_cg)
        % Conjugate Gradient (see below)
        t0 = tic;
        [u,iters] = cg(@matvec,f,tol,maxit);
        t1 = toc(t0);
    else
        % Matlab bicgstab
        t0 = tic;
        [u,flag,~,~,iters] = bicgstab(@matvec,f(:),tol,maxit);
        t1 = toc(t0);
        if (flag > 0)
            fprintf('Flag = %d\n',flag);
        end
        u = reshape(u,size(u0_shape));
    end
    
    % Compute errors
    e = abs(u_exact(2:N,2:N) - u(2:N,2:N));
    h2 = h^2;
    err(1,i) = sum(e(:))*h2/area;
    err(2,i) = sqrt(sum(e(:).^2)*h2/area);  
    err(3,i) = max(e(:));
        
    figure(1);
    clf;
    hd = surf(xm,ym,u);
    set(hd,'edgecolor','none');
    title('Solution');

    figure(2);
    loglog(iters,'.-','markersize',20);
    hold on;
    set(gca,'fontsize',18);
    itcount(i) = length(iters);
    title('Residual');
    
    tot_time(i) = t1;
    fprintf('%5d %10d %16.4e %16.4f\n',N,itcount(i),err(np,i),t1);
    
end
figure(3);
plot(nvec,itcount,'.-','markersize',20);
set(gca,'fontsize',18);
title('Iteration count');

figure(4);
loglog(nvec,err(np,:),'.','markersize',18);
set(gca,'fontsize',18);
s = polyfit(log(nvec),log(err(np,:)),1);
hold on;
loglog(nvec,exp(s(2))*nvec.^s(1),'-','linewidth',2);
str = sprintf('Conv. Rate %6.2f',s(1));
lh = legend(str);
title('Convergence rate');

figure(5);
loglog(nvec,tot_time,'.-','markersize',20);
set(gca,'fontsize',18);
title('Time (s)');

figure(6);
dofxits = (nvec(:).^2).*itcount(:);
y = dofxits./tot_time(:);
loglog(tot_time(:)./itcount(:),y,'.-','markersize',20);
set(gca,'fontsize',18);
title('Time (s)');
xlabel('Time','fontsize',16);
ylabel('DOFs x CG iterations/(seconds)','fontsize',16);



end

function [u,iters] = cg(matvec,f,tol,maxit)

global u0_shape;

u0 = zeros(size(u0_shape));

[mx,~] = size(u0);
N = mx-1;
iters = [];

r = f - matvec(u0);
u = u0;
p = r;
rn_old = mat_dot(r,r);
for k = 1:maxit
    Ap = matvec(p);
    alpha = rn_old/mat_dot(p,Ap);
    u = u + alpha*p;
    r = r - alpha*Ap;
    rn_new = mat_dot(r,r);
    iters(k) = rn_new;
    % fprintf('%4d %12.4e\n',k,rn_new);
    if (sqrt(rn_new) < tol)
        return;
    end
    b = rn_new/rn_old;
    p = r + b*p;    
    rn_old = rn_new;
end
fprintf('CG did not converge\n');

end

function Au = matvec(uin)

global b1 b2 u0_shape;

u = reshape(uin,size(u0_shape));

[mxp1,~] = size(u);
N = mxp1-1;
h = 1/N;

Au = zeros(N+1,N+1);

for i = 2:N
    for j = 2:N
        uxp = u(i+1,j) - u(i,j);
        uxm = u(i,j)   - u(i-1,j);
        uyp = u(i,j+1) - u(i,j);
        uym = u(i,j)   - u(i,j-1);
        Au(i,j) = (b1(i,j)*uxp - b1(i-1,j)*uxm + b2(i,j)*uyp ...
            - b2(i,j-1)*uym)/h^2;
    end
end

Au = reshape(Au,size(uin));

end

%{
function [b,bx,by] = beta_coeffs(x,y)

% Grady's Beta
b = 1 + x.*y;
bx = y;
by = x;

end
%}
    
function f = frhs(x,y)

y2 = y.^2;
y3 = y2.*y;
y4 = y3.*y;
y5 = y4.*y;

x2 = x.^2;
x3 = x2.*x;
x4 = x3.*x;
x5 = x4.*x;

f  = exp(x.*y).*(x5.*(y-1).*y2 ...
    - x4.*y.*(y2 - 7*y + 4) ...
    + x3.*(y5 - y4 - 6*y2 + 12*y - 3) ...
    - x2.*(y5 - 7*y4 + 6*y3 + 8*y - 5) ...
    - 2*x.*(2*y4 - 6*y3 + 4*y2 + 1) ...
    + y.*(-3*y2 + 5*y - 2));

end


function s = mat_dot(v1,v2)

s = sum(sum(v1.*v2));

end
