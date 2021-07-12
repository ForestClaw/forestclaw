function plots_paper()

global beta
global r0_torus
global xc0_torus yc0_torus zc0_torus
global alpha

pi2 = 2*pi;

[example,~,refine_pattern, alpha, beta, init_radius,~,~,tr,pr] = read_vars();

print_plots = true;

r0_torus = init_radius;


if example == 0 || example == 2
    set_init_point(1);
else
    set_init_point(3);
end
[xdata,ydata,zdata] = torus_init(512);

tfinal = 0.2;
save_data = true;
if (save_data)
    f = fopen('ring.dat','w');
    fprintf(f,'%d\n',length(xdata));
    fprintf(f,'%24.16f\n',tfinal);
    for i = 1:length(xdata)
        fprintf(f,'%24.16f %24.16f %24.16f\n',xdata(i),ydata(i),zdata(i));
    end
    fclose(f);
end    

% Plot path
plot_path = true;
if (plot_path)    
    hold on;
    tol = 0.02;

    % Plot paths
    if example == 0
        th = linspace(pi2*(0.25 - 1/32),pi2*(0.25+1/32),200);
        xpath = cos(th);
        ypath = sin(th);
        zpath = (alpha+0.01)*ones(size(xpath));
    elseif example == 1
        phi = pi2*linspace(0.125,0.375,200);
        r = 1 + (alpha+tol)*cos(phi);
        xpath = r.*cos(pi/2);
        ypath = r.*sin(pi/2);
        zpath = (alpha+tol)*sin(phi);
    end
    plot3(xpath,ypath,zpath,'k','linewidth',4);
   
    

    if example == 0 || example == 2
        set_init_point(1,alpha+tol);
    else
        set_init_point(3,alpha+tol);
    end
    fill_disk(xdata,ydata,zdata);   
    
    
    plot3(xc0_torus,yc0_torus,zc0_torus,'k.','markersize',30);
    
    if exist('ring.out','file')
        data = load('ring.out');
        xc = tr(1) + (tr(2)-tr(1))*data(:,1);
        yc = pr(1) + (pr(2)-pr(1))*data(:,2);
        [xth,yth,zth] = mapc2m_torus(xc,yc,alpha,beta);
                
        if example == 0 || example == 2
            set_init_point(2,alpha+tol);
        else
            set_init_point(4,alpha+tol);
        end
        % plot3(xth,yth,zth,'k','linewidth',3);
        fill_disk(xth,yth,zth);
    end        
    hold off;
end

v = [0, 65];
view(v);

showgridlines
showpatchborders
setpatchborderprops('linewidth',2)

setslicecolor('white')
o = findobj('type','colorbar');
if ishandle(o)
    delete(o);
end

daspect([1,1,1])
axis off
delete(get(gca,'title'))


if print_plots
    dir = '/Users/calhoun/projects/ForestClaw/papers/sisc_2014/figs_torus/tmp';
    if refine_pattern == 1
        print('-dpdf','-fillpage',[dir,'/torus_1.pdf']);
    else
        print('-dpdf','-fillpage',[dir,'/torus_2.pdf']);
    end
end


end

function plot_disk(p,N)

set_init_point(p);
[xdata,ydata,zdata] = torus_init(N);
hold on;
fill_disk(xdata,ydata,zdata);

end

function fill_disk(xdata,ydata,zdata)

global xc0_torus yc0_torus zc0_torus

h = [];
for i = 1:length(xdata)-1
    x = [xdata(i:i+1);xc0_torus];
    y = [ydata(i:i+1);yc0_torus];
    z = [zdata(i:i+1);zc0_torus];
    hold on;
    h(i) = fill3(x,y,z,[1,1,1]*0.8);
    set(h(i),'edgecolor','none','facealpha',0.7);
end    
    
plot3(xdata,ydata,zdata,'k','linewidth',2);

end

function set_init_point(disk_position,alpha1)

global xc0_torus yc0_torus zc0_torus beta alpha

if nargin < 2
    alpha1 = alpha;
end
   

% Center of the disk
if (disk_position == 1)
    th = 2*pi*(0.25 + 1.d0/32.d0);
    phi = pi/2;
elseif disk_position == 2
    th = 2*pi*(0.25 - 1.d0/32.d0);
    phi = pi/2;
elseif disk_position == 3
    th = pi/2;
    phi = pi/4;
elseif disk_position == 4
    th = pi/2;
    phi = 3*pi/4;
end
x0 = (1 + alpha1*cos(phi))*cos(th);
y0 = (1 + alpha1*cos(phi))*sin(th);
z0 = alpha1*sin(phi);

xc0_torus = x0;
yc0_torus = y0;
zc0_torus = z0;

end

% Get points that make up zero-level contour
function [xdata,ydata,zdata] = torus_init(N)

global beta alpha

xc = linspace(0,1,N);
yc = linspace(0,1,N);
[xcm,ycm] = meshgrid(xc,yc);

fm = torus_surface(xcm,ycm);

h = contourc(xc,yc,fm,[0,0]);

m = length(h);
xdata = nan(m,1);   % Final length < m
ydata = nan(m,1);
zdata = nan(m,1);

st_idx = 1;
k = 1;
while (1)
  % cval = h(1,st_idx);
  next_length = h(2,st_idx);

  xl = h(1,st_idx+1:st_idx+next_length);
  yl = h(2,st_idx+1:st_idx+next_length);

  [x0,y0,z0] = mapc2m_torus(xl(:),yl(:),alpha,beta);
  
  [x0,y0,z0] = fix_points(x0,y0,z0);
  
  xdata(k:k+next_length-1) = x0;
  ydata(k:k+next_length-1) = y0;
  zdata(k:k+next_length-1) = z0;
  k = k + next_length;
  st_idx = st_idx + next_length + 1;
  if (st_idx > length(h))
      break;
  end
end

% Remove extra nan's at the end of the array.
m = isnan(xdata);
xdata(m) = [];
ydata(m) = [];
zdata(m) = [];


end

function f = torus_surface(xc,yc)

global r0_torus xc0_torus yc0_torus zc0_torus beta alpha;

[xp,yp,zp] = mapc2m_torus(xc,yc,alpha,beta);

r2 = (xp-xc0_torus).^2 + (yp-yc0_torus).^2 + (zp - zc0_torus).^2 - r0_torus^2;
f = r2;

end

function [x0,y0,z0] = fix_points(X,Y,Z)

global alpha beta

pi2 = 2*pi;

% [~,~,~,~, ~,~,~,~,theta_range,phi_range] = read_vars();

th0 = atan2(Y,X);
m = ~isreal(th0);
if sum(m(:)) > 0
    fprintf('fix_points : th0 contains imaginary entries\n');
    keyboard;
end

m = th0 < 0;
th0(m) = th0(m) + pi2;

% tr1 = pi2*theta_range(1);
% tr2 = pi2*theta_range(2);
% xm = (th0-tr1)/(tr2-tr1);

phi = asin(Z/alpha);
m = ~isreal(phi);
if sum(m(:)) > 0
    fprintf('fix_points : phi contains imaginary entries\n');
    keyboard;
end

r = sqrt(X.^2 + Y.^2);

m1 = (r <= 1) & (phi > 0);
m2 = (r <= 1) & (phi <= 0);

phi(m1) =  pi - phi(m1);
phi(m2) = -pi - phi(m2);

m = phi < 0;
phi(m) = phi(m) + pi2;

r = 1 + alpha*cos(phi);
x0 = r.*cos(th0);
y0 = r.*sin(th0);
z0 = alpha*sin(phi);

end
