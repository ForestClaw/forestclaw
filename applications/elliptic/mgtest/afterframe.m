ax = 0;
bx = 1;
ay = 0;
by = 1;
s = 1e-2;
axis([ax-s bx+s ay-s by+s])
daspect([1 1 1]);
axis off;

showpatchborders;
setpatchborderprops('linewidth',1);

view(2);

rhs_choice = 3;
if rhs_choice == 3
    th = linspace(0,2*pi,500);
    r = 0.25;
    hold on;
    plot(r*cos(th)+0.5,r*sin(th)+0.5,'k','linewidth',2);
    hold off;
end

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'swirl0000','png');
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
clear mapc2m_fivepatch;
