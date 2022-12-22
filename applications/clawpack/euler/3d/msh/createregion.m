function createregion(R)
%
% CREATEREGION creates a topographic region of MSH for plotting or computing
%
%     CREATEREGION, by itself, allows the user to interactively select a
%     region of the Mt. St. Helens DEM file and save the region to a file.
%     This file can then be read for plotting or for constructing a grid
%     mapping for computational purposes.
%
%     CREATEREGION(R) uses the coarsening ratio R, which determines how
%     resolved the graphics will be.  On some machines, the rendering for
%     R = 1 may be quite slow. The default is R = 4.  This coarsening
%     factor only affects the visual display, not the actual grid
%     information that is stored in the user defined file.
%
%     This code assumes the Digital Elevation Map (DEM) file, named
%     "mshblalb.asc" and supplied by the USGS, is located in the
%     current directory.
%
%     See also READREGION, PLOTREGION, MAPC2P_MSH.
%

close all;

global msh_exit_val;

if (nargin == 0)
  R = 4;  % Refinement factor
end

% Read in data set from the file "mshblalb.asc"
ds_big = create_dataset();

% All relevant parameters are stored in the current figure window
% application
% space.
store_R(R);

% This plots  a top-down view of MSH, with no region selected.
plotregion(ds_big,R);

x = ds_big.x;
y = ds_big.y;
z = ds_big.msh_zhigh+100;

xlow = x(1);
xhigh = x(end);
ylow = y(1);
yhigh = y(end);

% Mcells and ncells are the number of cells in each direction in
% the full topographic map.
mcells = length(x) - 1;
ncells = length(y) - 1;

% this is to get around a matlab bug which returns the wrong
% axis position if the data aspect ratio has been set.
f = (yhigh - ylow)/(xhigh - xlow);
pa = [50 50 400 fix(f*400)];

% This function is recursively called and allows the user to
% draw and resize a region on the MSH data set.
dmove(xlow,xhigh,mcells,ylow,yhigh,ncells,z,pa)

if (msh_exit_val == 0)
  return;
end;

box = get_box();

ds_save = create_dataset(box);
fname = input('Input file name (''.region'' will be appended ) : ','s');
fname = [fname,'.region'];
fprintf('Storing region information in file ''%s''\n',fname);
write_dataset(fname,ds_save);

return;
% end of main routine CREATEREGION.


% -----------------------------------------------------------------
function view_data(ds,R)

% Internal graphics routine

% This routine plots a second view of the selected region.  This 3d
% view allows the user to visualize their selected region.

box = get_box();
ds = create_dataset(box);
R = get_R();

figure(2)
prect = get(gcf,'position');
set(gcf,'position',[20 prect(2:4)]);
plotregion(ds,R);
h = findobj(gca,'type','light');
if (isempty(h))
  camlight;
end;
view(3);
figure(1);

% -------------------------------------------------------
function store_R(R)
udata = get(gcf,'userdata');
udata.R = R;
set(gcf,'userdata',udata);

% -------------------------------------------------------
function R = get_R(R)
udata = get(gcf,'userdata');
R = udata.R;

% ---------------------------------------------------------
function write_dataset(fname,ds)

% Write the data set for use later

fid = fopen(fname,'w');

fprintf(fid,'%12d    Number of grid points in x\n',ds.mcells+1);
fprintf(fid,'%12d    Number of grid points in y\n',ds.ncells+1);

fprintf(fid,'%12.3f    xlow\n',ds.xmin);
fprintf(fid,'%12.3f    ylow\n',ds.ymin);

fprintf(fid,'%12.3f    xhigh\n',ds.xmax);
fprintf(fid,'%12.3f    yhigh\n',ds.ymax);

fprintf(fid,'%12.3f    Elevation low\n',ds.msh_zlow);
fprintf(fid,'%12.3f    Elevation high\n',ds.msh_zhigh);

fprintf(fid,'%12.3f     dx\n',ds.dx);
fprintf(fid,'%12.3f     dy\n',ds.dy);

fprintf(fid,'\n');
for i = 1:ds.mcells+1,
  for j = 1:ds.ncells+1,
    fprintf(fid,'%12.3f\n',ds.data(i,j));
  end;
  fprintf(fid,'\n');
end;
fclose(fid);


% --------------------------------------------------------
function datablock = create_dataset(box)

figure(1);
if (isappdata(gcf,'msh') == 0)
  % Read the data if it isn't already stored
  read_data;
end

% Get data from current figure window
msh = getappdata(gcf,'msh');

if (nargin == 0)
  ilo = 1;
  ihi = msh.mcells+1;
  jlo = 1;
  jhi = msh.ncells+1;
else
  if (isempty(box))
    datablock = [];
    return;
  end
  ilo = box.ilo;
  ihi = box.ihi;
  jlo = box.jlo;
  jhi = box.jhi;
end

ds = msh.data(ilo:ihi, jlo:jhi);
x = msh.x(ilo:ihi);
y = msh.y(jlo:jhi);

datablock.mcells = length(x)-1;
datablock.ncells = length(y)-1;
datablock.xmin = min(x);
datablock.ymin = min(y);
datablock.xmax = max(x);
datablock.ymax = max(y);
datablock.msh_zlow = msh.zmin;
datablock.msh_zhigh = msh.zmax;
datablock.dx = msh.cellsize;
datablock.dy = msh.cellsize;
datablock.data = ds;
datablock.x = x;
datablock.y = y;

% -------------------------------------------------------------
function read_data();

% Read raw data from DEM file and fix it up so we can use it
% Data will be stored in current figure window as app data.
% This data will be reference by 'create_dataset', above

varname = {'ncols','nrows','xllcorner','yllcorner', ...
    'cellsize','NODATA_value'};

fid = fopen('mshblalb.asc','r');

for i = 1:length(varname),
  [data,info] = fscanf(fid,'%s',1);
  if (strcmp(varname{i},char(data)) == 0)
    error('Variable name not read correctly');
  end;
  [data,info] = fscanf(fid,'%d',1);
  eval([varname{i},' = data;']);
  str = ['fprintf(''%s %d\n'',varname{i},',varname{i},');'];
  eval(str);
end;

fprintf('Reading data...\n');
[data,info] = fscanf(fid,'%f',nrows*ncols);
fprintf('Done reading data...\n');
data(data == NODATA_value) = nan;
msh_data = reshape(data,ncols,nrows);

x = xllcorner + (0:ncols)*cellsize;
y = yllcorner + (0:nrows)*cellsize;


% Create data set that is 32 blocks in the x direction by 23 blocks in the y
% direction. Each block is 40 grid values, and so represents 40*30 = 1200m
% of data.

% Choose xstart:xend, ystart:yend so that the data doesn't contain
% any NaN's (i.e. value -9999)
xstart = 7;
xend = 1286;
ystart = 8;
yend = 927;

% New data set without the -9999 values.
msh_new = msh_data(xstart:xend,ystart:yend);

% This gives us the x and y coordinates of the data.
x = x(xstart:xend);
y = y(ystart:yend);

msh.mcells = length(x)-1;
msh.ncells = length(y)-1;
msh.x = x;
msh.y = y;
msh.xmin = min(x);
msh.xmax = max(x);
msh.ymin = min(y);
msh.ymax = max(y);
msh.cellsize = cellsize;
msh.zmin = min(min(msh_new));
msh.zmax = max(max(msh_new));
msh.data = msh_new;

figure(1);
setappdata(gcf,'msh',msh);

% -------------------------------------------------------------
% Functions for creating the box
% This function is called recursively when various mouse
% events are initiated.
% -------------------------------------------------------------
function dmove(arg,xhigh,mcells,ylow,yhigh,ncells,z,pa)

global msh_exit_val;

if (nargin == 8)
  xlow = arg;
  arg = 'Initialize';
else
  box = get_box;
  [x,y,z,dx,dy,mcells,ncells] = get_info();
  xdir_rev = strcmp(get(gca,'xdir'),'reverse') == 1;
  ydir_rev = strcmp(get(gca,'ydir'),'reverse') == 1;
end;

switch arg
case 'Initialize'
   x = linspace(xlow,xhigh,mcells+1);
   y = linspace(ylow,yhigh,ncells+1);
   dx = (xhigh - xlow)/mcells;
   dy = (yhigh - ylow)/ncells;

   ha = findobj('type','axes');
   if (isempty(ha))
     error('You must first set up an image to crop');
   end;

   % Don't store anything until we have set up the axis.
   store_info(x,y,z,dx,dy);

   box = [];
   store_box(box);

   ab = [];
   store_ab(ab);

   s = 1;
   store_stepsize(s);

   % Hardwire position of axis to get around Matlab bug.
   axis_position_rect = pa;
   store_pa(axis_position_rect);
   set(ha,'units','pixels');
   set(ha,'position',axis_position_rect);

   clc;
   dmove input;

case 'input'

   clc;
   fprintf('CREATEREGION\n');
   fprintf('n        Create a new box\n');
   fprintf('f        Fix the aspect ratio of box\n');
   fprintf('s        Fix minimum stepsize\n');
   fprintf('m        Move/resize the selected box\n');
   fprintf('i        Information about the current box\n');
   fprintf('k        Enter Keyboard mode\n');
   fprintf('v        View current data set in 3d in new window\n');
   fprintf('r        Change display resolution\n');
   fprintf('w        Write out data file and quit\n');
   fprintf('q        Quit\n');
   fprintf('\n');

   ch = input('Input choice : ','s');

   figure(1);
   switch lower(ch)
   case 'n'
      dmove newbox;
   case 'f'
      dmove fix;
   case 's'
      dmove stepsize;
   case 'm'
      dmove moveresize;
   case 'i'
      dmove information;
   case 'k'
      keyboard;
      dmove input;
   case 'v'
      dmove view;
   case 'r'
      dmove resolution;
   case 'w'
      dmove write;
   case 'q'
      dmove quit;
   otherwise
      dmove input;
   end;

case 'newbox'

   str = ['-- Use the mouse to drag a rectangle over the desired ',...
	 'region.\n-- To move the rectangle, click the rectangle ',...
	 'with the left mouse button,\n   then click second time to new ',...
	 'location.\n-- Use yellow crosses to resize rectangle'];

   print_message(str,1);
   % store_ab([]);
   % store_stepsize(1);
   box = get_region(x,y,dx,dy,mcells,ncells);
   draw_box(box,x,y,z);
   dmove input;

case 'view'
   view_data;

   dmove input;

case 'resolution'
   R = input('Input new refinement factor for 3d view : ');
   if (~isempty(R))
     if (R > 0)
       store_R(R);
     end;
   end;
   dmove view;

case 'fix'
    fprintf('\n');
    fprintf('Input vector in the form [WIDTH, HEIGHT]. Example : [2,3]\n');
    fprintf('\n');
    ab = input('Input new aspect ratio : ');
    while (length(ab) ~= 2)
      fprintf('\n');
      fprintf('Input aspect ratio as a vector ''[WIDTH, HEIGHT]'' \n');
      ab = input('Input new aspect ratio : ');
    end;
    store_ab(ab);

    if (~isempty(box))
      fixed_pt = 3;  % Fix upper left corner and adjust lower right
      box = fix_aspect_ratio(box,mcells,ncells,fixed_pt);
      draw_box(box,x,y,z);
      store_box(box);
    end;

    dmove input;

case 'stepsize'
    s = input('Input stepsize : ');
    store_stepsize(s);

    if (~isempty(box))
      fixed_pt = 3;
      box  = fix_aspect_ratio(box,mcells,ncells,fixed_pt);
      draw_box(box,x,y,z);
      store_box(box);
    end;

    dmove input;

case 'information'
    ds_new = create_dataset(box);
    if (isempty(ds_new))
      str = 'No box is currently selected or created.';
      print_message(str,0);
    else
      fprintf('\n\n');
      fprintf('Number of cells in x direction     : %5i\n',ds_new.mcells);
      fprintf('Number of cells in y direction     : %5i\n',ds_new.ncells);
      fprintf('Longitude range                    : [%6.0f, %6.0f]\n',...
	  ds_new.xmin, ds_new.xmax);
      fprintf('Latitude range                     : [%6.0f, %6.1f]\n',...
	  ds_new.ymin, ds_new.ymax);
      fprintf('Size of region (lat x long, (km))  : %4.1f x %4.1f\n',...
	  (ds_new.xmax - ds_new.xmin)/1000, (ds_new.ymax - ...
	  ds_new.ymin)/1000);
      fprintf('Height minimum (meters)            : %4.1f\n',...
	  min(min(ds_new.data)));
      fprintf('Height maximum (meters)            : %4.1f\n',...
	  max(max(ds_new.data)));
      fprintf('Aspect ratio                       : %6.2f\n',...
	  ds_new.mcells/ds_new.ncells);
      fprintf('\n');
      print_message('',2);
    end;

    dmove input;

case 'write'
  msh_exit_val = 1;
  return;

case 'quit'

  msh_exit_val = 0;
  return;

case 'anchor1'
   irect = get_box_rect(box,x,y);

   % Upper right point is fixed point
   fpt = [irect(1) + irect(3); irect(2) + irect(4)];
   rbbox(irect,fpt);
   point2 = get(gca,'currentpoint');

   p1 = get_fixed_point(1,box,x,y,xdir_rev,ydir_rev);
   p2 = point2(1,1:2);

   box = make_new_box(p1,p2,x,y,dx,dy,mcells,ncells);

   draw_box(box,x,y,z);
   store_box(box);

   dmove nothing;

case 'anchor2'
   irect = get_box_rect(box,x,y);

   % lower right point is fixed point.
   fpt = [irect(1) + irect(3); irect(2)];

   rbbox(irect,fpt);
   point2 = get(gca,'currentpoint');

   p1 = get_fixed_point(2,box,x,y,xdir_rev,ydir_rev);
   p2 = point2(1,1:2);

   box = make_new_box(p1,p2,x,y,dx,dy,mcells,ncells);

   draw_box(box,x,y,z);
   store_box(box);

   dmove nothing;

case 'anchor3'
   irect = get_box_rect(box,x,y);

   % Lower left point is fixed point
   fpt = [irect(1) irect(2)];

   rbbox(irect,fpt);
   point2 = get(gca,'currentpoint');

   p1 = get_fixed_point(3,box,x,y,xdir_rev,ydir_rev);
   p2 = point2(1,1:2);

   box = make_new_box(p1,p2,x,y,dx,dy,mcells,ncells);

   draw_box(box,x,y,z);
   store_box(box);

   dmove nothing;

case 'anchor4'
   irect = get_box_rect(box,x,y);

   % Upper left point is the fixed point
   fpt = [irect(1) irect(2)+irect(3)];

   rbbox(irect,fpt);
   point2 = get(gca,'currentpoint');

   p1 = get_fixed_point(4,box,x,y,xdir_rev,ydir_rev);
   p2 = point2(1,1:2);

   box = make_new_box(p1,p2,x,y,dx,dy,mcells,ncells);

   draw_box(box,x,y,z);
   store_box(box);

   dmove nothing;

case 'bdown'
   point1 = get(gca,'currentpoint');
   k = waitforbuttonpress;
   point2 = get(gca,'currentpoint');

   p1 = point1(1,1:2);
   p2 = point2(1,1:2);

   % Get distance to move the box
   dp = p2 - p1;
   idist = fix(dp./[dx dy]);

   irange = [1-box.ilo mcells + 1 - box.ihi];
   jrange = [1-box.jlo ncells + 1 - box.jhi];

   ip = max([min([idist(1),irange(2)]),irange(1)]);
   jp = max([min([idist(2),jrange(2)]),jrange(1)]);

   box.ilo = box.ilo + ip;
   box.ihi = box.ilo + box.m;
   box.jlo = box.jlo + jp;
   box.jhi = box.jlo + box.n;

   draw_box(box,x,y,z);
   store_box(box);

   dmove nothing;

case 'nothing'
   % ----------------------------------
end;


% -----------------------------------------------------
function box = get_region(x,y,dx,dy,mcells,ncells)

delete_square;
delete_anchors;
delete_grid;
delete_markers;


k = waitforbuttonpress;

% Starting point for box
point1 = get(gca,'CurrentPoint');

% Rubberband box
rbbox;

% Ending point for box
point2 = get(gca,'CurrentPoint');

p1 = point1(1,1:2);
p2 = point2(1,1:2);

box = make_new_box(p1,p2,x,y,dx,dy,mcells,ncells);

% -----------------------------------------------------
function box = make_new_box(p1,p2,x,y,dx,dy,mcells,ncells)


% Now get indices in array
x1 = p1(1);
y1 = p1(2);
x2 = p2(1);
y2 = p2(2);

% sort points to get lower left and upper right
llx = min([x1 x2]);
lly = min([y1 y2]);
urx = max([x1 x2]);
ury = max([y1 y2]);

% Get indices for box
[cv1,box.ilo] = min(abs(llx-x));
[cv2,box.jlo] = min(abs(lly-y));
[cv3,box.ihi] = min(abs(urx-x));
[cv4,box.jhi] = min(abs(ury-y));

box.m = box.ihi - box.ilo;
box.n = box.jhi - box.jlo;

s = sign(p2 - p1);
if (s(1) > 0 & s(2) > 0)
  fixed_pt = 1;
elseif (s(1) > 0 & s(2) <= 0)
  fixed_pt = 2;
elseif (s(1) <= 0 & s(2) <= 0)
  fixed_pt  = 3;
elseif (s(1) <= 0 & s(2) > 0)
  fixed_pt = 4;
end;

box = fix_aspect_ratio(box,mcells,ncells,fixed_pt);

store_box(box);

% -----------------------------------------------------------------
function draw_box(box,x,y,z)

delete_square;
delete_anchors;
delete_grid;
delete_markers;

% Note that these must go in this order
draw_grid(box,x,y,z)
draw_square(box,x,y,z)
draw_markers(box,x,y,z);
draw_anchors(box,x,y,z);

% -----------------------------------------------------------------
function draw_square(box,x,y,z)

xlow  = x(box.ilo);
xhigh = x(box.ihi);

ylow  = y(box.jlo);
yhigh = y(box.jhi);

hr = line('xdata',[xlow xlow  xhigh xhigh xlow],...
          'ydata',[ylow yhigh yhigh ylow  ylow],...
	  'zdata',[z z z z z]);
set(hr,'color','r','linewidth',2);
set(hr,'tag','square');
set(hr,'buttondownfcn','dmove bdown');

% -----------------------------------------------------------------
function delete_square()

hr = findobj('tag','square');
if (~isempty(hr))
  delete(hr);
end;

% -----------------------------------------------------------------
function draw_anchors(box,x,y,z)

if (strcmp(get(gca,'xdir'),'reverse') == 1)
  xlow = x(box.ihi);
  xhigh = x(box.ilo);
else
  xlow = x(box.ilo);
  xhigh = x(box.ihi);
end;

if (strcmp(get(gca,'ydir'),'reverse') == 1)
  ylow = y(box.jhi);
  yhigh = y(box.jlo);
else
  ylow = y(box.jlo);
  yhigh = y(box.jhi);
end;

anchor(1) = line('xdata',xlow, 'ydata',ylow, 'buttondownfcn','dmove anchor1');
anchor(2) = line('xdata',xlow, 'ydata',yhigh,'buttondownfcn','dmove anchor2');
anchor(3) = line('xdata',xhigh,'ydata',yhigh,'buttondownfcn','dmove anchor3');
anchor(4) = line('xdata',xhigh,'ydata',ylow, 'buttondownfcn','dmove anchor4');
set(anchor,'zdata',z,'marker','+','color','y','linewidth',3, ...
    'markersize',10);
set(anchor,'tag','anchors');

% -----------------------------------------------------------------
function delete_anchors()
anchors = findobj('tag','anchors');
if (~isempty(anchors))
  delete(anchors);
end;

% -----------------------------------------------------------------
function draw_grid(box,x,y,z)

ab = get_ab();
if (isempty(ab))
  return;
end;

a = ab(1);
b = ab(2);

if (a/b == 1)
  % for aspect ratio of 1, draw a 2x2 grid
  str = 'Aspect ratio reduced to [2,2]';
  print_message(str,0);
  if (fix(box.m/2) ~= box.m/2)
    mt = (box.m-1)/2;
    nt = (box.n-1)/2;
  else
    mt = box.m/2;
    nt = box.n/2;
  end;
else
  mt = box.m/a;
  nt = box.n/b;
end;

xv = x(box.ilo:mt:box.ihi);
yv = y(box.jlo:nt:box.jhi);

for i = 1:length(xv),
  grid(i) = line('xdata',[xv(i) xv(i)],'ydata',[yv(1) yv(end)],...
      'zdata',[z z]);
end;

for j = 1:length(yv),
  grid(length(xv) + j) = line('xdata',[xv(1) xv(end)],...
      'ydata',[yv(j) yv(j)],'zdata',[z z]);
end;

cl = [1 1 1]*0.8;
set(grid,'linestyle','-','color',cl,'tag','grid','linewidth',1);

% ----------------------------------------------------------------
function delete_grid()
grid = findobj('tag','grid');
if (~isempty(grid))
  delete(grid);
end;

% -----------------------------------------------------------------
function draw_markers(box,x,y,z)

step = get_stepsize();

if (step < 5)
  % this just doesn't look good
  % return;
end

xv = x(box.ilo:step:box.ihi);
yv = y(box.jlo:step:box.jhi);

for i = 2:length(xv)-1,
  markers_x(i-1) = line('xdata',[xv(i) xv(i)],'ydata',[yv(1) yv(end)],...
      'zdata',[z z]);
end;

for j = 2:length(yv)-1,
  markers_y(j-1) = line('xdata',[xv(1) xv(end)],...
      'ydata',[yv(j) yv(j)],'zdata',[z z]);
end;

cl = [1 1 1]*0.8;
if (step <= 10)
  msize = 1;
elseif (step <= 30)
  msize = 3;
else
  msize = 5;
end;

% Note that we set the buttondown function so that if the user
% hits a marker  instead of the actual box, the box can still
% be moved.

set([markers_x markers_y],'linestyle','none','color',cl,...
    'tag','markers','marker','+','markersize',msize,...
    'buttondownfcn','dmove bdown');

if (step == 1)
  set([markers_x markers_y],'visible','off');
end

% ----------------------------------------------------------------
function delete_markers()
markers = findobj('tag','markers');
if (~isempty(markers))
  delete(markers);
end;

% ----------------------------------------------------------
function box = fix_aspect_ratio(box,mcells,ncells,fixed_pt)

switch fixed_pt
case 1
    s = [1 1];
    istart = box.ilo;
    jstart = box.jlo;
case 2
    s = [1 -1];
    istart = box.ilo;
    jstart = box.jhi;
case 3
    s = [-1 -1];
    istart = box.ihi;
    jstart = box.jhi;
case 4
    s = [-1 1];
    istart = box.ihi;
    jstart = box.jlo;
end;


ab = get_ab;
step = get_stepsize;

if (isempty(ab))
  m_new = floor(box.m/step)*step;
  n_new = floor(box.n/step)*step;

  iend = istart + s(1)*m_new;
  jend = jstart + s(2)*n_new;
else

  % a,b should be relatively prime
  [a,b] = reduce(ab(1),ab(2));

  store_ab([a b]);

  p = 0;
  while (step*a*p <= box.m)
    p = p + 1;
  end;

  q = 0;
  while (step*b*q <= box.n)
    q = q + 1;
  end;

  r = min([p q]) - 1;

  m_mdpt = mcells/2 + 1;
  n_mdpt = ncells/2 + 1;

  while (true)
    iend = istart + s(1)*a*r*step;
    jend = jstart + s(2)*b*r*step;
    done = (abs(iend - m_mdpt) <= mcells/2 + 1 & ...
	abs(jend - n_mdpt) <= ncells/2 + 1);
    if (done)
      break;
    end;
    r = r - 1;
    if (r == 0)
      error('could not find any pair that works!');
    end;
  end;
end;

switch fixed_pt
case 1
    box.ihi = iend;
    box.jhi = jend;
case 2
    box.ihi = iend;
    box.jlo = jend;
case 3
    box.ilo = iend;
    box.jlo = jend;
case 4
    box.ilo = iend;
    box.jhi = jend;
end;

box.m = box.ihi - box.ilo;
box.n = box.jhi - box.jlo;

% -----------------------------------------------------------------
function [a,b] = reduce(a,b)

if (a/b == 1)
  a = 1;
  b = 1;
  return;
end;

while (true)
  af = factor(a);
  bf = factor(b);

  common_factors = intersect(af,bf);
  if (isempty(common_factors))
    break;
  end
  for i = 1:length(common_factors)
    a = a/common_factors(i);
    b = b/common_factors(i);
  end;
end;


% -----------------------------------------------------------------
function irect = get_box_rect(box,x,y)

% Get corners of the box in axis (pixel) units
% Units are in pixels,relative to the lower left corner of
% the figure window.  note that we need to use 'pa',the
% axis position rectangle, that we explicitly set in the
% initialization phase.
% We use this rectangle to resize the box later.

pa = get_pa();

p1 = [x(box.ilo) y(box.jlo)];
p2 = [x(box.ihi) y(box.jhi)];

% Get location of box in figure coordinates (pixels).
if (strcmp(get(gca,'xdir'),'reverse') == 1)
  p1f(1) = fix(pa(1) + ((x(end) - p2(1))/(x(end) - x(1)))*pa(3));
  p2f(1) = fix(pa(1) + ((x(end) - p1(1))/(x(end) - x(1)))*pa(3));
else
  p1f(1) = fix(pa(1) + ((p1(1) - x(1))/(x(end) - x(1)))*pa(3));
  p2f(1) = fix(pa(1) + ((p2(1) - x(1))/(x(end) - x(1)))*pa(3));
end;

if (strcmp(get(gca,'ydir'),'reverse') == 1)
  p1f(2) = fix(pa(2) + ((y(end) - p2(2))/(y(end) - y(1)))*pa(4));
  p2f(2) = fix(pa(2) + ((y(end) - p1(2))/(y(end) - y(1)))*pa(4));
else
  p1f(2) = fix(pa(2) + ((p1(2) - y(1))/(y(end) - y(1)))*pa(4));
  p2f(2) = fix(pa(2) + ((p2(2) - y(1))/(y(end) - y(1)))*pa(4));
end;


% Set up rectangle in figure units.
d = p2f - p1f;
ilof = p1f(1);
jlof = p1f(2);
irect = [ilof jlof abs(d)];

% ----------------------------------------------------------------
function p1 = get_fixed_point(num,box,x,y,xdir_rev,ydir_rev)

switch num
case {1,2}
  if (xdir_rev)
    p1(1) = x(box.ilo);
  else
    p1(1) = x(box.ihi);
  end;
case {3,4}
  if (xdir_rev)
    p1(1) = x(box.ihi);
  else
    p1(1) = x(box.ilo);
  end;
end;

switch num
case {1,4}
  if (ydir_rev)
    p1(2) = y(box.jlo);
  else
    p1(2) = y(box.jhi);
  end;
case {2,3}
  if (ydir_rev)
    p1(2) = y(box.jhi);
  else
    p1(2) = y(box.jlo);
  end;
end;



% ----------------------------------------------------------------
% Functions for getting and storing information
% ----------------------------------------------------------------

% -----------------------------------------------------------------
function store_box(box)

udata = get(gcf,'Userdata');
udata.box = box;
set(gcf,'userdata',udata);

% -------------------------------------------------------------------
function box = get_box()

udata = get(gcf,'userdata');
box = udata.box;

% -----------------------------------------------------------------
function store_ab(ab)

udata = get(gcf,'Userdata');
udata.ab = ab;
set(gcf,'userdata',udata);

% -------------------------------------------------------------------
function ab = get_ab()

udata = get(gcf,'userdata');
ab = udata.ab;

% -----------------------------------------------------------------
function store_stepsize(s)

udata = get(gcf,'Userdata');
udata.stepsize = s;
set(gcf,'userdata',udata);

% -------------------------------------------------------------------
function s = get_stepsize()

udata = get(gcf,'userdata');
s = udata.stepsize;

% -----------------------------------------------------------------
function store_pa(pa)

udata = get(gcf,'Userdata');
udata.pa = pa;
set(gcf,'userdata',udata);

% -------------------------------------------------------------------
function pa = get_pa()

udata = get(gcf,'userdata');
pa = udata.pa;

% --------------------------------------------------------
function store_info(x,y,z,dx,dy)

info.x = x;
info.y = y;
info.z = z;
info.dx = dx;
info.dy = dy;
udata = get(gcf,'UserData');
udata.info = info;
set(gcf,'userdata',udata);

% --------------------------------------------------------
function [x,y,z,dx,dy,mcells,ncells] = get_info()

udata = get(gcf,'UserData');
info = udata.info;

x = info.x;
y = info.y;
z = info.z;
dx = info.dx;
dy = info.dy;
mcells = length(x) - 1;
ncells = length(y) - 1;


% -----------------------------------------------------------
function print_message(str,flag)

if (flag == 0)
  fprintf('\n');
  fprintf(sprintf('%s\n\n',str));
  input('Hit <enter> to continue...');
elseif(flag == 1)
  fprintf('\n');
  fprintf(sprintf('%s\n\n',str));
else
  fprintf('\n\n');
  input('Hit <enter> to continue...');
end;
