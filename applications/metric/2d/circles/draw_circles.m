function [circle_h, diag_h, box_h] = draw_circles(box)


radius = box.r1;

th = linspace(0,2*pi,201);
xdata = radius.*cos(th);
ydata = radius.*sin(th);

xlow = box.xlow;
xhi = box.xhi;
ylow = box.ylow;
yhi = box.yhi;

xcircle = (xdata+1)*(xhi-xlow)/2 + xlow;
ycircle = (ydata+1)*(yhi-ylow)/2 + ylow;

% Add circle
circle_h = line('xdata',xcircle,'ydata',ycircle,'zdata',0*xdata+0.0,'color','r','linewidth',3);

% Add diagonals
diag_h(1) = line('xdata',[xlow,xhi],'ydata',[ylow,yhi],'zdata',[0,0],'color','b','linewidth',2);
diag_h(2) = line('xdata',[xlow,xhi],'ydata',[yhi,ylow],'zdata',[0,0],'color','b','linewidth',2);

% Add Box
box_h = line('xdata',[xlow xhi xhi xlow xlow],'ydata',[ylow ylow yhi yhi ylow],...
    'zdata',[0 0 0 0 0],'color','k','linewidth',2);
