function parallel_colors(p,x,y,z,q)

qcolors = q;

% Use logarithmic scale for errors
m = q > 0;
qlog = q;
qlog(m) = log10(q(m));

npmax = 4;  % total number of procs
n = 64;  % length of jet colormap;
mpirank = getmpirank();


% Map q linearly into colormap; constrain entire color map into
% range [qmin, qmax]. 
qmin = -100;
qmax = -8;

qm = qlog(m);    % log of error, e.g -20, -10, -4
cidx = zeros(size(qm));

% Everything less than qmin will take processor color
m1 = qm < qmin;
cidx(m1) = mpirank + 1;

% Everything between qmin and qmax, map linearly into 'jet'
m2 = qmin <= qm & qm <= qmax;
if (~isempty(qm))
    cidx(m2) = floor((n-1)/(qmax-qmin)*(qm(m2) - qmin)) + 1;  
end
cidx(m2) = cidx(m2) + npmax;

m3 = qm > qmax;
cidx(m3) = npmax + n + 1;

qcolors = q;
qcolors(m) = cidx;
qcolors(~m) = mpirank + 1;

set(p,'cdata',qcolors);
set(p,'cdatamapping','direct');
set(p,'facecolor','flat');



end