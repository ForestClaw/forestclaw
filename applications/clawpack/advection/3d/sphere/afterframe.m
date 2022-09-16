setviews;

yrbcolormap;
caxis([0 1]);

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

showpatchborders;
setpatchborderprops('linewidth',1);
% showgridlines;

parms = read_vars();

R = parms.maxelev;
axis([-1-R,1+R,-1-R,1+R,-1-R,1+R]);
set(gca,'box','on')
view(3)

% camlight;


NoQuery = 0;
prt = false;
if (prt)
    filename = sprintf('sphere_%0.4d.png',Frame);
    disp(filename);
    print('-dpng',filename);
end

shg