function plot_tetondam()

close all;

% plot_topo('topos/TetonDamSmallHiRes_orig.topo','g');
% plot_topo('topos/TetonDamLargeLowRes.topo','b');
% plot_topo('topos/TetonDamFloodPlain.topo','r');
[ax,bx,ay,by] = plot_topo('topos/TetonLarge.topo','w');

% Lower left   ( -112.34626736,  43.18013542)
% Upper right  ( -111.26428819,  43.95986458)

axis([-112.34626736,-111.26428819,43.18013542,43.95986458]);
shg

end

function [ax,bx,ay,by] = plot_topo(fname,c)
% PLOT_TOPO plots the topo from file FNAME


% Top of the cinder cone.  Pass these values into plot_feature
% to get the height of the feature.
% xp = 3.0206e+03;
% yp = 1.1689e+04;

% behind the dam
% xp =  2.6732e+03;
% yp = 1.1942e+04;

% c = 'w';  % Color of the topo
hold on;
[p,ax,bx,ay,by] = plot_feature(fname,c);
hold on;

%fprintf('Height at input location : %12.4f\n',hp);

daspect([1,1,1e4]);

% axis([ax bx ay by]);

camlight;
setviews;
view(vtop);
shg

end


function [p,ax,bx,ay,by,hpout] = plot_feature(fname,c,xp,yp)

fid = fopen(fname);

ncols = fscanf(fid,'%d',1); fscanf(fid,'%s',1);
nrows = fscanf(fid,'%d',1); fscanf(fid,'%s',1);
xll = fscanf(fid,'%g',1);   fscanf(fid,'%s',1);
yll = fscanf(fid,'%g',1);   fscanf(fid,'%s',1);
dx = fscanf(fid,'%g',1);    fscanf(fid,'%s',1);
nodata = fscanf(fid,'%g',1); fscanf(fid,'%s',1);
T = fscanf(fid,'%g',nrows*ncols);
fclose(fid);

% --------------------------------
ax = xll;
ay = yll;

bx = ax + dx*(ncols-1);
by = ay + dx*(nrows-1);

x = linspace(ax,bx,ncols);
y = linspace(ay,by,nrows);

T = reshape(T,ncols,nrows);
T = fliplr(T);
% T = reshape(T,nrows,ncols)';
% T = fliplr(T);


[xm,ym] = meshgrid(x,y);

nodata = find(T == nodata);
T(nodata) = nan;
c = T;
c(~nodata) = 0;
c(nodata) = 1;

p = surf(xm,ym,T');
set(p,'cdata',c');
set(p,'edgecolor','none');

colormap([1 0 0; 1 1 1]);
set(p,'cdatamapping','direct');

plot_water = true;
if (plot_water)
    hold on;
    
    pc = path_coords();

    x0 = -111.5391666666667;
    x1 = xll + (ncols-1)*dx;
    h0 = 1580;
    h1 = 1720;
    f = @(x,y) h0 + (h1-h0)/(x1-x0)*(x-x0);
    xw = pc(:,1);
    yw = pc(:,2);
    zw = f(xw,yw);
    h = fill3(xw,yw,zw,'b');
    set(h,'facealpha',0.7);
    
    % water = patch(surf2patch(xw,yw,zw));
    % set(water,'facecolor','b');
    % set(water,'edgecolor','k','linewidth',2);
    % set(water,'facealpha',0.7);
end
plot_zero = false;
if (plot_zero)
    hold on;
    % Add plane
    xll = 22155.63;
    yll  = 10531.91;
    y0 = 11689;
    h0 = 100;
    h1 = 230;
    h(ym < y0) = nan;
    set(gcf,'renderer','opengl');
    xw = [0 5400; 0 5400];
    yw = [0 0; 25000 25000];
    zw = [0 0; 0 0];


    water = patch(surf2patch(yw,xw,zw));
    set(water,'facecolor','b');
    set(water,'edgecolor','k','linewidth',2);
    set(water,'facealpha',0.7);
end


fprintf('Min height  : %12.4f\n',min(T(:)));
fprintf('Max height  : %12.4f\n',max(T(:)));

zlim([min(T(:)), max(T(:))]);


if (nargin > 2)
    % find height of given location (xp,yp)
    hp = interp2(xm,ym,T,xp,yp,'linear');
    if (nargout > 5)
        hpout = hp;
    end
end


end

function pc = path_coords()

pc= [ ...
-111.5329602010579,43.90768251711403,0 
-111.5425694320295,43.91291520283087,0 
-111.5265517148,43.92985173935236,0 
-111.4838301614854,43.93342747699161,0 
-111.4526923402719,43.92275343741503,0 
-111.4229733394132,43.92102950267466,0 
-111.3888345969115,43.942256354582464,0 
-111.3547747442403,43.96012412188881,0 
-111.3290981366292,43.94067321304157,0 
-111.3055308323582,43.92361145688078,0 
-111.2787647823537,43.93643579897589,0 
-111.2483464534973,43.94374241316081,0 
-111.2512048905199,43.93181718339741,0 
-111.3004354725815,43.91434621684393,0 
-111.3310923535663,43.92063116805807,0 
-111.3574777595737,43.93748955926343,0 
-111.3734430072187,43.92390465892466,0 
-111.3830291449872,43.91639601995375,0 
-111.4070584125505,43.91053159674234,0 
-111.427329896465,43.90372138216379,0 
-111.4497771036826,43.90847035302123,0 
-111.4470054376565,43.8931939107954,0 
-111.4340292718339,43.88128609194279,0 
-111.4328310870467,43.87453080137639,0 
-111.452551818065,43.87284083784596,0 
-111.4655146366021,43.88648791247213,0 
-111.474254800272,43.89924049270922,0 
-111.4743486255008,43.91523470848247,0 
-111.4878354873365,43.92113549058154,0 
-111.5039146384135,43.91948128966648,0 
-111.5182133156775,43.91450172601051,0]; 

pc(end+1,:) = pc(1,:);
end
