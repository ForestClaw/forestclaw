function daspect_brick()
clf


ax = -112.36171859;
bx = -111.25911794;
ay = 43.59190493;
by = 43.97790751;

% Meter's
d1 = compute_distances([ax bx],[ay,ay]);
d2 = compute_distances([ax,ax],[ay,by]);

fprintf('%20s %12.4f\n','Horizontal distance',d1/1e3);
fprintf('%20s %12.4f\n','Vertical distance',d2/1e3);
fprintf('\n\n');

use_distance = true;
if (use_distance)
    ax = 0;
    bx = d1;
    ay = 0;
    by = d2;
end


plot([ax bx bx ax ax],[ay ay by by ay],'k','linewidth',2);
hold on;
daspect([1 1 1]);

% Aspect ratio : 2.85645927 (actual) vs. 2.8421052631578947 (approx)
mi = 2;
mj = 1;

xe = linspace(ax,bx,mi+1);
ye = linspace(ay,by,mj+1);


for i = 1:(mi+1)
    plot([xe(i) xe(i)],[ay by],'k');
end

for j = 1:(mj+1)
    plot([ax,bx],[ye(j) ye(j)],'k');
end

set(gca,'fontsize',16);
xlabel('longitude','fontsize',16);
ylabel('latitude','fontsize',16);

% GeoClaw finegrid resolution : 
R = [1,2,4,4,4];
r = 1;
mx = 64;
my = 32;
dx = (bx-ax)/mx;
dy = (by-ay)/my;
fprintf('%7s %16s %16s %8s %8s\n','level','dx','dy','mx','my');
fprintf('%s\n',char(double('-')*ones(1,59)));
for i = 1:length(R)    
    dx = dx/R(i);
    dy = dy/R(i);
    mx = mx*R(i);
    my = my*R(i);
    if (use_distance)
        fprintf('%7d %16.2f %16.2f %8d %8d\n',i,dx,dy,mx,my);
    else
        fprintf('%7d %16.4e %16.4e %8d %8d\n',i,dx,dy,mx,my);
    end
end
fprintf('%s\n',char(double('-')*ones(1,59)));
fprintf('\n');

% ForestClaw
minlevel = 0;
maxlevel = 10;
mx0 = 16;
mx = mi*mx0*2^minlevel;
my = mj*mx0*2^minlevel;
dx = (bx-ax)/mx;
dy = (by-ay)/my;
fprintf('%7s %16s %16s %8s %8s\n','level','dx','dy','mx','my');
fprintf('%s\n',char(double('-')*ones(1,59)));
for l = minlevel:maxlevel    
    if (use_distance)
        fprintf('%7d %16.2f %16.2f %8d %8d\n',l,dx,dy,mx,my);
    else
        fprintf('%7d %16.4e %16.4e %8d %8d\n',l,dx,dy,mx,my);
    end
    dx = dx/2;
    dy = dy/2;
    mx = mx*2;
    my = my*2;
end
fprintf('%s\n',char(double('-')*ones(1,59)));

shg

end

function d = compute_distances(lon,lat)

% lon = [ax bx bx ax ax];
% lat = [ay ay by by ay];

fprintf('\n');
for i = 1:1,
    lat1 = lat(i);
    lat2 = lat(i+1);
    lon1 = lon(i);
    lon2 = lon(i+1);
    
    dlon = lon2-lon1;
    dlat = lat1-lat2;
    
    R = 6367.5e3;
    a = (sind(dlat/2))^2 + cosd(lat1) * cosd(lat2) * (sind(dlon/2))^2 ;
    c = 2 * atan2( sqrt(a), sqrt(1-a) ) ;
    d = R * c;
%     fprintf('%10s %d : %8.0f\n','Distance',i,d);
end


end
