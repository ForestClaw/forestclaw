function convert2latlong(topo_num)

if topo_num == 1
    % Small Hi Res topo
    ncols = 834;
    nrows = 180;
    ncells = 30;

    xll = 22155.63;
    yll = 10531.91;
    
    ax = -111.684952;
    bx = -111.373222;

    ay = 43.889109;
    by = 43.937331;

    dx = (bx-ax)/(ncols-1);
    dy = (by-ay)/(nrows-1);
    
    fprintf('dx = %12.6e\n',dx);
    fprintf('dy = %12.6e\n',dy);

    fprintf('Topo distances for topo 1\n');
    compute_distances(ax,bx,ay,by);
elseif (topo_num == 2)
    % Use topo 1 data to get coordinates for big topo
    
    xll = 22155.63;  % Units in flood plane coordinates
    th_ll = -111.684952;
    
    xlr = 47145.63;
    th_lr = -111.373222;
    
    slope = (th_ll - th_lr)/(xll-xlr);
    th_xll = slope*(0-xll) + th_ll;
    th_xlr = slope*(48000 - xll) + th_ll;
    
    yll = 10531.91;
    th_ll = 43.889109;
    
    yur = 10531.91 + 30*179;
    th_ur = 43.937331;

    slope = (th_ll - th_ur)/(yll-yur);
    th_yll = slope*(0-yll) + th_ll;
    th_yur = slope*(17450 - yll) + th_ll;

    fprintf('Topo coordinates for Flood plain topo\n');
    
    fprintf('lower left  : %16.8f %16.8f\n',th_xll,th_yll);
    fprintf('upper right : %16.8f %16.8f\n',th_xlr,th_yur);

   
elseif (topo_num == 3)  
    % TetonDamLarge (in Cartesian coordinates)
    ncols = 3031;
    nrows = 1431;
    ncells = 30;

    xll = 387718.130435160303;
    yll = 4826352.473530278541;
    
    % from gdalinfo -json
    ll_coords = [
        [
        -112.4 
        43.98
        ] 
        [
        -112.3907344,
        43.5816825
        ] 
        [
        -111.2298244
        43.5899112
        ]
        [
        -111.2313561
        43.9883432
        ]
        [
        -112.4 
        43.98
        ]];

    ll_coords = reshape(ll_coords,2,5);
    
    long_lim = ll_coords(1,:);
    lat_lim = ll_coords(2,:);
    
    fprintf('Topo distances for topo 3\n');
    compute_distances(long_lim,lat_lim);
    
elseif (topo_num == 4)
    % Large topo (in LatLong coordinates)
   
    
    ncols = 4180;
    nrows = 1464;
    ncells = 0.000277729665;

    
    ll_coords =  [
        [
          -112.3907344 
          43.9883432
        ] 
        [
          -112.3907344
          43.5816825
        ] 
        [
          -111.2298244
          43.9883432
        ]
        [
          -111.2298244
          43.5816825
        ]
        [
          -112.3907344
          43.9883432
        ]];
    
    ll_coords = reshape(ll_coords,2,5);
    
    long_lim = ll_coords(1,[1 2 4 3 5]);
    lat_lim = ll_coords(2,[1 2 4 3 5]);
    
    fprintf('Topo distances for topo 4\n');
    compute_distances(long_lim,lat_lim);

    
end

% % Gauge data
% ay = 43.913661;
% by = 43.916382;
% 
% ax = -111.623926;
% bx = -111.620150;
% 
% fprintf('Plant distances\n');
% compute_distances(ax,bx,ay,by);




end



function compute_distances(lon,lat)

% lon = [ax bx bx ax ax];
% lat = [ay ay by by ay];

fprintf('\n');
for i = 1:4,
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
    fprintf('%10s %d : %8.0f\n','Distance',i,d);
end

end

