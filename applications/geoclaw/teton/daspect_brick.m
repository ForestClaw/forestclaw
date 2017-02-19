clf


ax = -112.36171859;
bx = -111.25911794;
ay = 43.59190493;
by = 43.97790751;

% Aspect ratio : 2.85645927 (actual) vs. 2.8421052631578947 (approx)
mi = 17;
mj = 6;

plot([ax bx bx ax ax],[ay ay by by ay],'k','linewidth',2);
hold on;
daspect([1 1 1]);

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
mx = 54;
dx = (bx-ax)/mx;
fprintf('%7s %16s\n','level','Res.');
fprintf('%s\n',char(double('-')*ones(1,24)));
for i = 1:length(R)    
    dx = dx/R(i);
    fprintf('%7d %16.4e\n',i,dx);
end
fprintf('%s\n',char(double('-')*ones(1,24)));
fprintf('\n');

% ForestClaw
minlevel = 0;
maxlevel = 6;
mx = mi*8*2^minlevel;
my = mj*8*2^minlevel;
dx = (bx-ax)/mx;
dy = (by-ay)/my;
fprintf('%7s %16s %16s\n','level','dx','dy');
fprintf('%s\n',char(double('-')*ones(1,41)));
for l = minlevel:maxlevel
    fprintf('%7d %16.4e %16.4e\n',l,dx,dy);
    dx = dx/2;
    dy = dy/2;
end
fprintf('%s\n',char(double('-')*ones(1,41)));

shg
