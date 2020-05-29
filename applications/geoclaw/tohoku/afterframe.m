% Patches
hidegridlines(1:6);
showpatchborders(1:7);
% hidepatchborders(5);
setpatchborderprops('linewidth',2);
% set(gca,'zlim',[-20 1]);   % Need so that all patchborders show up

colormap(parula);
colorbar;
tol = -0.8;
c1 = -0.1;
c2 = 0.1;
caxis([c1,c2]);

%{
cv1 = linspace(qmin,c1,11);
cv2 = linspace(c2,qmax,11);
cv = [cv1,cv2];
drawcontourlines(cv);
%}

set(gca,'color',[101,67,33]/255)

hold on;
% add_gauges();
add_regions(t);
hold off;

fprintf('%20s %12.4e\n','qmin',qmin);
fprintf('%20s %12.4e\n','qmax',qmax);


% Axes
axis([132 210 9 53])
daspect([1 1 1]);
set(gca,'fontsize',16);
axis([174.7508,  209.4175,   13.2162,  32.7718]);

title(sprintf('Tohoku : t = %.2f',t),'fontsize',18);

NoQuery = 0;
prt = false;
MaxFrames = 1000;
if (prt)
    filename = sprintf('bowl%04d.png',Frame);
    fprintf('Print file %s\n',filename);
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
