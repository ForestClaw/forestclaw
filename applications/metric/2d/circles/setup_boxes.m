function [M,n] = setup_boxes()

% -----------------------------------
% Set up boxes
% -----------------------------------

create_boxes();

fid = fopen('boxes.dat','r');
data = textscan(fid,'%f',2,'CommentStyle','%');
nbox = data{1}(1);
M = data{1}(2);
% [nbox,M] = deal(data{1}(1),data{1}(2),data{1}(3));

% data = textscan(fid,'%f',1,'CommentStyle','%');

for i = 1:nbox,
  data = textscan(fid,'%f %f %f %f',1,'CommentStyle','%');
  for m = 1:2,
    smallend{i}(m) = data{m}(1);
  end
  box_size(i) = data{3}(1);
  box_radius(i) = data{4}(1);
end
fclose(fid);

for i = 1:nbox,
  n = box_size(i);
  r = box_radius(i);
  bigend = smallend{i} + [n,n]-1;
  add_box(M,smallend{i},bigend,r,i,n);
end;

function create_boxes()

global boxes;
boxes = struct('xlow',-1,'xhi',1,'ylow',-1,'yhi',1,'r1',0.5,'size',2);

% --------------------------------------------------------
% Construct the boxes
% --------------------------------------------------------
function add_box(M,smallend,bigend,r1,i,n)

% CONSTRUCT_BOX sets up structure needed to include circle
%
% BOX = CONSTRUCT_BOX(M,SMALLEND,BIGEND,R) Sets up a box defined on an
% [MxM] grid, with lower left corner in cell SMALLEND and upper right corner
% in cell BIGEND.   The size of M only needs to be large enough to resolve
% your desired box layout.
%
% The circle enclosed by this box is specified by indicating a radius,
% 0 < R < 1.  The circle will then occupy (pi*R^2) fraction of the box.
%
% Remark : Do not make R too close to 1.  You should always have a buffer of
% a few grid cells between the enclosing box and the boundary of the circle.
%

global boxes;

h = 2/M;
% xgrid = linspace(-1,1,M+1);
% ygrid = linspace(-1,1,M+1);
xlow = -1 + (smallend(1) - 1)*h;
xhi = -1 + bigend(1)*h;
ylow = -1 + (smallend(2) - 1)*h;
yhi = -1 + bigend(2)*h;

% boxes(i) = ...
%     struct('smallend',smallend,'bigend',bigend,'xgrid',xgrid,'ygrid',ygrid,...
%     'xlow',xlow,'xhi',xhi,'ylow',ylow,'yhi',yhi,'h',h,'r1',r1);

boxes(i) = ...
    struct('xlow',xlow,'xhi',xhi,'ylow',ylow,'yhi',yhi,'r1',r1,'size',n);
