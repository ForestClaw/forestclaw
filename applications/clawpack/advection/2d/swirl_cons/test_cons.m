function [gridsout,solnout] = test_cons()

% ----------------------
% Physical parameters
% ----------------------
ax = 0;
bx = 1;

uvel = @(x) cos(2*pi*x) + 2;   % Velocity field

uinit = @(x) ones(size(x));    

% Time stepping
N = 5;
dt = 2.5e-3;
T = dt*N;

prt_soln = false; 

% ----------------------
% Numerical parameters 
% ----------------------
mx = 8;

% grid setup
sx = 1/mx;

se = linspace(0,1,mx+1);
sc = se(1:end-1) + sx/2;

levels = [3 3 3 4 4 4 4 3 3 3];
mg = length(levels);

grids = cell(1,mg);
soln = cell(N+1,1);

soln{1}.q = cell(1,mg);
xlower = ax;
sum0 = 0;
for m = 1:mg
    l = levels(m);
    w = (bx-ax)*2^(-l);

    dx = sx*w;    
    
    g.level = l;
    g.xlower = xlower;
    g.xupper = xlower + w;
        
    g.xe = linspace(g.xlower-dx,g.xupper+dx,mx+3);
    g.xc = g.xe(1:end-1) + dx/2;
    g.dx = dx;
    g.aux = uvel(g.xc);
    grids{m} = g;

    % Initial conditions
    q0 = uinit(g.xc);
    sum0 = sum0 + sum(q0(2:end-1))*g.dx;
    soln{1}.q{m} = q0;
    
    xlower = g.xupper;
end

fprintf('sum0 = %24.16f\n',sum0);

% Update
for n = 1:N
    solnq = fill_bc(grids,soln{n}.q);    
    rp = riemann(grids,solnq);
    soln{n+1}.q = cell(1,mg);
    soln{n+1}.qold = solnq;
    sumc_before = 0;
    for m = 1:mg
        g = grids{m};
        dx = g.dx;
        qn = soln{n}.q{m};  % No BCs
        qnp1 = qn(2:end-1) - dt/dx*(rp{m}.amdq(2:end) + rp{m}.apdq(1:end-1));
        sumc_before = sumc_before + sum(qnp1)*dx;
        
        soln{n+1}.q{m} = [nan qnp1 nan];   % Fill in BCs later
    end 
    [soln{n+1},delta,fdiff,rdiff] = apply_correction(grids,rp,dt,soln{n+1});
    sumc_after = 0;
    for m = 1:mg
        dx = grids{m}.dx;
        sumc_after = sumc_after + sum(soln{n+1}.q{m}(2:end-1))*dx;
    end
    
    fprintf('Frame %d : dt = %12.6f     t = %12.6f\n',n,dt,n*dt);
    fprintf('sum[0] =  %24.16f %24.16e\n',sumc_before, abs(sumc_before-sum0));
    fprintf('sum[0] =  %24.16f %24.16e\n',sumc_after, abs(sumc_after-sum0));
    fprintf('\n');
    fprintf('Delta  =  %24.16f %24.16f\n',delta(1), delta(2));
    fprintf('fdiff  =  %24.16f %24.16f\n',fdiff(1), fdiff(2));
    fprintf('rdiff  =  %24.16f %24.16f\n',rdiff(1), rdiff(2));
    fprintf('\n');
end

if (prt_soln)
    close all;
    plot_soln(grids,soln);
end

if (nargout)
    gridsout = grids;
    solnout = soln;
end
end

function [soln,delta,fdiff,rdiff] = apply_correction(grids, rp, dt, soln_new)

mg = length(grids);

soln = soln_new;

% Right edge of grid 1
ic = 3;

gc = grids{ic};
gf = grids{ic+1};

ql = soln_new.qold{ic}(end-1);
qr = soln_new.qold{ic+1}(1);   % ghost cell value
ul = gc.aux(end-1);
ur = gf.aux(1);

fmc = rp{ic}.amdq(end);
fmf = rp{ic+1}.amdq(1);
fdiff(1) = fmf-fmc;
rdiff(1) = ur*qr - ul*ql;
delta(1) = fdiff(1) + rdiff(1);
soln.q{ic}(end-1) = soln.q{ic}(end-1) - dt/gc.dx*delta(1);


% Left edge of grid mg
ic = 8;

gc = grids{ic};
gf = grids{ic-1};

ql = soln_new.qold{ic-1}(end);    % ghost cell value
qr = soln_new.qold{ic}(2);   
ul = gf.aux(end);
ur = gc.aux(2);

fpc = rp{ic}.apdq(1);
fpf = rp{ic-1}.apdq(end);
fdiff(2) = fpf - fpc;
rdiff(2) = ur*qr - ul*ql;
delta(2) = fdiff(2) + rdiff(2);
soln.q{ic}(2) = soln.q{ic}(2) - dt/gc.dx*delta(2);

end


function rp = riemann(grids,solnq)

mg = length(grids);

rp = cell(1,mg);
for m = 1:mg
    g = grids{m};
    q = solnq{m};
    aux = g.aux;
    
    ql = q(1:end-1);
    qr = q(2:end);
    ul = aux(1:end-1);
    ur = aux(2:end);
    
    r.s = ur;
    r.waves = (ur.*qr - ul.*ql)./ur;
    r.amdq = zeros(size(ur));
    r.apdq = ur.*qr - ul.*ql;    
    rp{m} = r;
end



end

function solnq = fill_bc(grids,solnq)

mg = length(grids);


for m = 1:mg
    g = grids{m};
    if m == 1
        gl = grids{end};
        gr = grids{m+1};
        ql = solnq{end};
        qr = solnq{m+1};
    elseif m == mg
        gl = grids{m-1};
        gr = grids{1};
        ql = solnq{m-1};
        qr = solnq{1};
    else
        gl = grids{m-1};
        gr = grids{m+1};
        ql = solnq{m-1};
        qr = solnq{m+1};
    end
    
    
    % Left grid
    if (g.level > gl.level)
        dxc = gl.dx;
        % left grid is coarser grid
        qc(1) = ql(end-2);
        qc(2) = ql(end-1);
        qc(3) = mean(solnq{m}(2:3));   % Fill in average grid value.
        slope = (qc(3) - qc(1))/(2*dxc);
        solnq{m}(1) = qc(2) + slope*dxc/4;        
    elseif (g.level == gl.level)
        % left grid is same level
        solnq{m}(1) = ql(end-1);
    else
        % left grid is finer grid
        solnq{m}(1) = mean(ql(end-2:end-1));
    end
    
    % right grid
    if (g.level > gr.level)
        % right grid is coarser grid
        dxc = gr.dx;
        qc(1) = mean(solnq{m}(end-2:end-1));  % Fill in average grid value.
        qc(2) = qr(2);
        qc(3) = qr(3);
        slope = (qc(3) - qc(1))/(2*dxc);
        solnq{m}(end) = qc(2) - slope*dxc/4;
    elseif (g.level == gr.level)
        % right grid is same level
        solnq{m}(end) = qr(2);
    else
        % right grid is finer grid
        solnq{m}(end) = mean(qr(2:3));
    end
        
    
    
    
    
end


end


function plot_soln(grids, soln)

N = length(soln);
mg = length(grids);

for n = 1:N
    clf;
    solnq = soln{n}.q;
    for m = 1:mg
        g = grids{m};
        plot(g.xc(2:end-1),solnq{m}(2:end-1),'ro');
        hold on;
    end
    plot([0 1],[0 0],'k-');
    axis([0 1 -1.1 3.1]);
    set(gca,'fontsize',16);
    shg;
    pause;    
end


end






