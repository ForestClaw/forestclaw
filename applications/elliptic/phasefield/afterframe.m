domain_length = 1.5;
%{
if (t < 0.1)
    dl = 1;
elseif (0.1 <= t && t < 0.5)
    dl = 2;
elseif (0.5 <= t && t < 1.5)
    dl = 3;
elseif (1.5 <= t && t < 1.75)    
    dl = 3.5;
elseif (1.75 <= t && t < 2.0)    
    dl = 4;
elseif (2.0 <= t && t < 2.5)    
    dl = 5.0;
elseif (2.5 <= t && t < 3)    
    dl = 5.5;
else 
    dl = domain_length/2;
end
%}
dl = domain_length/2;
userAxis = [-1,1,-1,1]*dl;
axis(userAxis);
grid off;
daspect([1 1 1]);
axis square;
colorbar;

if (mq == 1)
    caxis([-1,0]);
    % caxis([-1,-0.99])
elseif (mq == 2)
    caxis([0,1]);
    % drawcontourlines([0.5,0.5]);
end

%{
if (mq == 1)
    [xcm,ycm] = meshgrid(xcenter,ycenter);
    phi = reshape(amrdata(1).data(2,:),mx,my);
    hold on;
    [~,h_phi] = contour(xcm,ycm,phi,[0.1,0.5,0.9],'k');

    cv = linspace(-1,0,11);
    cv([1 end]) = [];
    [~,h_u(1)] = contour(xcm,ycm,q,cv,'k');        
    [~,h_u(2)] = contour(xcm,ycm,q,[-0.99,-0.99],'k--');
    [~,h_u(3)] = contour(xcm,ycm,q,[-0.999,-0.999],'k--');
    hold off;
    set(h_phi,'linewidth',2);
end
%}
    

fprintf('%10s %16.8f\n','qmin',qmin);
fprintf('%10s %16.8f\n','qmax',qmax);
fprintf('%10s %16.8e\n','qmax-qmin',qmax-qmin);

showpatchborders;

set(gca,'fontsize',16);

% create files for a movie;  Use 'convert2gif' to convert .jpg files
% to .gif files;  Use 'convert -delay 10 dendrite???.gif dendrite.gif'
% to create movie dendrite.gif.
NoQuery = 0;
MaxFrames = 14;
mmovie = mod(Frame,1) == 0;
if (mmovie)
    prtname = sprintf('dendrite%03d',Frame);
    print('-dpng',prtname);
end

clear afterframe;
