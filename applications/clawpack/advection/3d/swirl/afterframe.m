yrbcolormap
axis([0 1 0 1 0 1])
daspect([1 1 1]);


fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);


showpatchborders;
setpatchborderprops('linewidth',1);
setpatchbordercolor('k');

% showcubes;
parms = read_vars();
showslices;
if parms.example > 0
    hideslices('x');
    hideslices('y');
end    

showsurfs();

caxis([0,1]);

set(gca,'box','on');

%axis off;
axis on;

shg;