yrbcolormap
axis([0 1 0 1 0 1])
daspect([1 1 1]);

hideslices;
showslices('x',6);
showslices('y',6);
showslices('z',6);

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);


showpatchborders;
setpatchbordercolor('k');
% setpatchborderprops(1,'linewidth',2,'color','k');  % new version only
% setpatchborderprops(2,'linewidth',2,'color','k');  % new version only
% setpatchborderprops(3,'linewidth',2,'color','k');  % new version only

showcubes;
setcubecolor('r',1);
setcubecolor('b',2);
setcubecolor('k',3);
hidecubes(1:2);

caxis([0,1]);

showgridlines(1:2);

% cv = linspace(0,1,11);
% cv([1 end]) = [];
% drawcontourlines(cv);

h = surflight;

axis off;

shg;

clear afterframe;
