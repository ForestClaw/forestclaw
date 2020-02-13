function [r,p] = map1d(xgrid,ygrid,qgrid)
%
%  MAP1D is a user-defined function to create line plots from 2d/3d data
%
%      MAP1D is a user-defined function called by PLOTCLAW2 or PLOTCLAW3 to
%      convert 2d or 3d data into 1d data for a line plot.
%
%      For 2d data, the MAP1D function has the form
%
%           [r,q] = map1d(xgrid,ygrid,qgrid)
%
%      where XGRID, YGRID, and QGRID are MxN arrays of grid data, and R and
%      Q are 1d column vectors of user-defined length.
%
%      For 3d data, the MAP1D function has the form
%
%           [r,q] = map1d(xgrid,ygrid,zgrid, qgrid)
%
%      where XGRID, YGRID, ZGRID and QGRID are MxNx{ arrays of grid data,
%      R and Q are 1d column vectors of user-defined length.
%
%      This function will be called from PLOTCLAW2 or PLOTCLAW3 if
%      'UserMap1d' = 1.   This plotting parameter can be set in the file
%      SETPLOT2 or SETPLOT3.
%
%      Example :
%
%             % To plot q as a function of x only ;
%             function [r,q] = map1d(xgrid,ygrid,qgrid)
%
%             [m,n] = size(xgrid);
%             r = reshape(xgrid,m*n,1);
%             q = reshape(qgrid,m*n,1);
%
%      Example :
%
%           % Plot data interpolated to grid diagonal on [0,1]x[0,1]x[0,1]
%           % grid
%           function [r,q] = map1d(xgrid,ygrid,zgrid,qgrid)
%
%           r = linspace(0,0.1,1)*sqrt(3);
%           q = interp3(xgrid,ygrid,zgrid,qgrid,r);
%
%
%      NOTE : This help file can also be found by typing 'map1d_help' at
%      the Matlab prompt
%
%
%  See also SETPLOTSTYLE, PLOTFRAME1EZ, GETLEGENDINFO, SETPLOT.

[mx,my] = size(qgrid);

% [xp,yp,~] = mapc2m(xgrid,ygrid);
r = reshape(xgrid,mx*my,1);
p = reshape(qgrid,mx*my,1);
