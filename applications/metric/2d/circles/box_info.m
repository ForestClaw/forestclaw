function box_info(box,dx)

xlow = box.xlow;
xhi = box.xhi;
ylow = box.ylow;
yhi = box.yhi;
r1 = box.r1;

ndx = (xhi - xlow)/dx;
ndy = (yhi - ylow)/dx;

ts = 1e5;
ndx = round(ndx*ts)/ts;
ndy = round(ndy*ts)/ts;

rx = r1*ndx/2;
ry = r1*ndy/2;
rx = round(rx*ts)/ts;
ry = round(ry*ts)/ts;

fprintf('\n');
fprintf('%-40s %5g\n','Cells in the x direction',ndx);
fprintf('%-40s %5g\n','Cells in the y direction',ndy);
fprintf('%-40s %5g\n','Cells in half circle in x direction',rx);
fprintf('%-40s %5g\n','Cells in half circle in y direction',ry);
fprintf('\n');
