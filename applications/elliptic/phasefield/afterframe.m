domain_length = 12;
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
    caxis([-1.1,0]);
    % caxis([-1,-0.99])
elseif (mq == 2)
    caxis([0,1]);
    drawcontourlines([0.5,0.5]);
end
colormap(parula);

%%{
if (mq == 1)
    cv = linspace(0,1,11);
    cv([1 end]) = [];
    tau = 1e-2;
    cv = [tau/2, 0.5, 1-tau/2];
    cv2 = [tau,1-tau];
    de = linspace(0,1,mx+1);    
    dc = de(1:end-1) + 0.5/mx;
    pmin = 10000;
    pmax = 0;
    clear h_phi h_u;
    for k = 1:length(amrdata)
        xlow = amrdata(k).xlow;
        ylow = amrdata(k).ylow;
        dx = amrdata(k).dx;
        dy = amrdata(k).dy;
        xhi = xlow + mx*dx;
        yhi = ylow + my*dy;
        xc = xlow + dc*(xhi - xlow);
        yc = ylow + dc*(yhi - ylow);
        [xcm,ycm] = meshgrid(xc,yc);
        
        % Plot phi
        phi = reshape(amrdata(k).data(2,:),mx,my)';
        pmin = min([pmin;phi(:)]);
        pmax = max([pmax;phi(:)]);
        hold on;
        if (min(phi(:)) < max(phi(:)))
            [~,h_phi(1)] = contour(xcm,ycm,phi,[0.5, 0.5],'k','linewidth',2);
            hold on;
            % [~,h_phi(2)] = contour(xcm,ycm,phi,[tau/2, 1-tau/2],'k');
        end
        % m = ishandle(h_phi);
        % set(h_phi(m),'linewidth',2);

        %{
        % Plot u
        cv_u = linspace(-1,0,11);
        cv_u([1 end]) = [];
        q = reshape(amrdata(k).data(1,:),mx,my)';
        if (min(q(:)) < max(q(:)))
            [~,h_u(1,k)] = contour(xcm,ycm,q,cv_u,'k');
            [~,h_u(2,k)] = contour(xcm,ycm,q,[-0.99,-0.99],'k--');
            [~,h_u(3,k)] = contour(xcm,ycm,q,[-0.999,-0.999],'k--');
        end
        %}
    end
           
    hold off;
else
    cv = linspace(0,1,11);
    cv([1 end]) = [];
    drawcontourlines(cv);
end
%}
    

fprintf('%10s %16.8e\n','pmin',pmin);
fprintf('%10s %16.8e\n','pmax',pmax);
fprintf('%10s %16.8e\n\n','pmax-pmin',qmax-qmin);
fprintf('%10s %16.8f\n','umin',qmin);
fprintf('%10s %16.8f\n','umax',qmax);
fprintf('%10s %16.8e\n','qmax-qmin',qmax-qmin);

showpatchborders;
setpatchborderprops('linewidth',1);
% hidepatchborders(8);

set(gca,'fontsize',16);

% create files for a movie;  Use 'convert2gif' to convert .jpg files
% to .gif files;  Use 'convert -delay 10 dendrite???.gif dendrite.gif'
% to create movie dendrite.gif.
NoQuery = 1;
MaxFrames = 32;
mmovie = mod(Frame,1) == 0;
if (mmovie)
    if (length(amrdata) == 1)
        prefix = 'single';
    else
        prefix = 'amr';
        hidepatchborders;        
    end
    setcontourlineprops('linewidth',1);
    fs = [16,16];
    dpi = 256;
    set(gcf,'paperunits','inches');
    set(gcf,'papersize',fs);
    set(gca,'position',[0.08 0.08 0.84 0.84]);
    set(gcf,'paperposition',[0 0 fs]);
    set(gcf,'GraphicsSmoothing','off')
    resstr = sprintf('-r%d',dpi);
    fname = sprintf('dendrite_%s_%02d.png',prefix,Frame);
    fprintf('Printing grid ''%s''\n',fname);
    print('-dpng',resstr,fname);
end

clear afterframe;
