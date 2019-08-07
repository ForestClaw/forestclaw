function plot_refine_contours(mx,maxlevel,t,u,v,refine_threshold)

% Define function and parameters
r0 = 0.2;
Hsmooth = @(r) (tanh(r/0.02d0) + 1)/2.d0;
qinit = @(r) Hsmooth(r + r0) - Hsmooth(r - r0);

x0 = 0.5;
y0 = 0.5;

% Clean up previous contour lines
h = findobj('tag','refine_region');
if ~isempty(h)
    delete(h);
end


% Draw contours
N = mx*2^maxlevel;

s = linspace(0,1,N+1);

x = mod(s - u*t,1);
y = mod(s - v*t,1);
[xm,ym] = meshgrid(x,y);

rm = sqrt((xm-x0).^2 + (ym-y0).^2);

q0 = qinit(rm);

rf = refine_threshold;
[~,h(1)] = contour(s,s,q0,[rf,rf],'k','linewidth',2);
hold on;

[~,h(2)] = contour(s,s,q0,[1-rf,1-rf],'k','linewidth',2);
set(h,'tag','refine_region');

daspect([1,1,1])

shg;

end


