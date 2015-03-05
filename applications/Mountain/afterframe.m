yrbcolormap;
showpatchborders(1:10);
caxis([0,1])

view(2);
axis([0 5000 0 2000]);
axis on;
daspect([1 1 1]);
set(gca,'fontsize',16,'fontweight','bold');

delete(get(gca,'title'));
delete(get(gca,'xlabel'));
delete(get(gca,'ylabel'));

mnt = true;
if (mnt)
    m = load('mountain.dat');
    xm = linspace(min(m(:,1)),max(m(:,1)),2000);
    ym = pchip(m(:,1),m(:,2),xm);    
    xm = [0; xm(:); xm(end); 0];
    ym = [0; ym(:); 0; 0];
    hold on;
    hmount_fill = fill(xm,ym,'k');
    set(hmount_fill,'visible','on','facecolor',[0 100/256 0]);
    hmount_line = line('xdata',xm,'ydata',ym,'linewidth',3);
    setpatchborderprops('linewidth',1);
    % colormap([0.6 0.7 0.9]);
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
