axis image;
daspect([1 1 1]);
axis off;

plot_choices = {'partition','solution','mesh'};
plot_type = plot_choices{2};
yrbcolormap;
qlo = 0.1;
qhi = 1.1;
under_label = sprintf('%3.1f - %7.1e',qlo,qlo-qmin);
over_label = sprintf('%3.1f + %7.1e',qhi,qmax-qhi);
fprintf('%6s %12s\n','qmin',under_label);
fprintf('%6s %12s\n\n','qmax',over_label);

if (ShowUnderOverShoots)
    colorbar_underover(under_label,over_label);
end;

% showgridlines(1:4);
showpatchborders(1:8);
setpatchborderprops(1:8,'linewidth',1);
view([-51.5,0]);

setviews;
% set(gcf,'visible','off');

NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 50;
  filename = framename(Frame,'mesh000','png');
  fprintf('Printing file %s\n',filename);
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
clear mapc2m_pillow;
