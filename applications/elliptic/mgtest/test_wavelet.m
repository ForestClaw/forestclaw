function test_wavelet(L)

global N0 lmax;

if (nargin == 0)
    L = 8;
end

N0 = 1;
N = N0*2^L;    

lmax = L;

h = 2/N;
xe = linspace(0,2,N+1);
xc = xe(1:end-1) + h/2;

% sin(2.*pi*x) + 0.4*sin(15*pi*x)*max(sin(2.*pi*x), 0);
f = @(x) sin(2*pi*x) + 0.4*sin(15*pi*x).*max([sin(2*pi*x);zeros(size(x))]);
% f = @(x) Hsmooth(x-0.5) - Hsmooth(x-1.5);
%f = @(x) abs(x-1) < 0.5;
% f = @(x) exp(-60*(x-1).^2);
% f = @(x) sin(8*pi*x).*exp(-20*(x-1).^2);



% Get fine grid solution (for plotting)
xf = linspace(0,2,2001);
qf = f(xf);


% Compute average function (pointwise would work too)
q = compute_average(f,xc);
% q = f(xc);

% Compute tranform
w = wavelet_transform(q);

close all

s = 2;
for l = 0:L
    nlevel = 2^l*N0;
    hlevel = 2/nlevel;
    xelevel = linspace(0,2,nlevel + 1);
    xclevel = xelevel(1:end-1) + hlevel/2;
    
    figure(1)
    % Plot average q (exact average)
    plot(xf,0*xf - s*l,'k');
    hold on;    
    plot(xf,f(xf)-s*l,'linewidth',2);
    plot(xclevel,w(l+1).q - s*l,'k.-','markersize',17);
    
    figure(2);
    % Plot transform
    plot(xf,0*xf - s*l,'k');
    hold on;    
    plot(xf,f(xf)-s*l,'linewidth',2);
    plot(xclevel,w(l+1).wt - s*l, 'k.-','markersize',17);    
end
yticks = -(L:-1:0)*s;
figure(1)
set(gca,'ytick',yticks);
set(gca,'yticklabels',L:-1:0);
title('Average solution down','fontsize',16);

figure(2)
set(gca,'ytick',yticks);
set(gca,'yticklabels',L:-1:0);
title('Wavelet transform','fontsize',16);

% Plot inverse transform
figure(3);
N = 2^L;
h = 2/N;
xe = linspace(0,2,N + 1);
xc = xelevel(1:end-1) + h/2;

q2 = inverse_wavelet(w);
plot(xf,qf,'r','linewidth',2);
hold on;
plot(xc,q2,'k.-','markersize',17);

% figure(4)
% semilogy(xc,abs(q2-q),'.-','markersize',17);
fprintf('Diff = %12.4e\n',norm(q2-q,Inf));


shg


end

function q = compute_average(f,xc)

% Coarsen the function
N = length(xc);
h = 2/N;

Nf = 200;
w = h/2;
hf = h/Nf;

for i = 1:N
    xi = xc(i);
    xef = linspace(xi-w,xi+w,Nf+1);
    xcf = xef(1:end-1) + hf/2;
    
    qf = f(xcf);
    q(i) = sum(qf)/Nf;
end

end

function w = wavelet_transform(q)

global lmax;

for l = lmax:-1:1
    qr = restriction(q);
    qp = prolongation(qr);
    
    w(l+1).q = q;
    w(l+1).qp = qp;
    w(l+1).wt = q-qp;
    q = qr;
end

w(1).q = qr;
w(1).qp = qp;
w(1).wt = 0;


end

function q = inverse_wavelet(w)

global lmax;

qcoarse = w(1).q;
for l = 1:lmax
    qp = prolongation(qcoarse);
    q = w(l+1).wt + qp; 
    qcoarse = q;
end

end


function qr = restriction(q)

qr = (q(1:2:end) + q(2:2:end))/2;

end

function qp = prolongation(q)

% Reinterpolate value
N = length(q);
h = 2/N;

qp = zeros(1,2*N);
xp = [-h, 0, h];

% Assume periodic boundary conditions
q_ext = [q(end), q, q(1)];
for i = 1:N
    yp = q_ext(i:i+2);
    p = polyfit(xp,yp,2);
    qp((i-1)*2 + 1) = polyval(p,-h/4);
    qp((i-1)*2 + 2) = polyval(p,+h/4);
end

end

function H = Hsmooth(r)

H = (tanh(r/0.01) + 1)/2;

end