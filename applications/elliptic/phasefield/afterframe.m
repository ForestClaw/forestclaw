domain_length = 12;
%{
if (t < 0.1)
    dl = 1;
elseif (0.1 <= t && t < 0.5)
    dl = 2;
elseif (0.5 <= t && t < 1.0)
    dl = 2.5;
elseif (1.0 <= t && t < 1.5)
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
dl = domain_length;
userAxis = [0,1,0,1]*dl;
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
    tau = 5e-3;
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
            [~,h_phi(k)] = contour(xcm,ycm,phi,[0.5 0.5],'k','linewidth',2);
            hold on;
            % Plot the transpose to check for symmetry.
            contour(ycm,xcm,phi,[0.5 0.5],'r','linewidth',2);
            % [~,h_phi(2)] = contour(xcm,ycm,phi,[tau/2, 1-tau/2],'k');
        end
        % m = ishandle(h_phi);
        % set(h_phi(m),'linewidth',2);

        plot_u = true;
        if (plot_u)
            cv_u = linspace(-1,0,11);
            cv_u([1 end]) = [];
            q = reshape(amrdata(k).data(1,:),mx,my)';
            if (min(q(:)) < max(q(:)))
                [~,h_u(1,k)] = contour(xcm,ycm,q,cv_u,'k');
                % [~,h_u(2,k)] = contour(xcm,ycm,q,[-0.99,-0.99],'k--');
                tol = 1e-3;    % Plot extent of thermal field.
                [~,h_u(3,k)] = contour(xcm,ycm,q,[-1+tol,-1+tol],'k-');
                tol = 1e-6;    % Plot extent of thermal field.
                [~,h_u(3,k)] = contour(xcm,ycm,q,[-1+tol,-1+tol],'k-');
            end
        end
    end
           
    hold off;
else
    cv = linspace(0,1,11);
    cv([1 end]) = [];
    drawcontourlines(cv);
end
%}

if (Frame <= 32)
    plot_fishpack = true;
else
    plot_fishpack = false;
end
if (plot_fishpack)
    % ax = copyobj(gca,gcf);
    hold on;
    [amr_fp,t] = readamrdata_forestclaw(2,Frame,'./fishpack/');
    mx_fp = amr_fp.mx;
    my_fp = amr_fp.my;
    xlow = amr_fp.xlow;
    ylow = amr_fp.ylow;
    dx = amr_fp.dx;
    dy = amr_fp.dy;
    xhi = xlow + mx_fp*dx;
    yhi = ylow + my_fp*dy;
    xe = linspace(xlow,xhi,mx_fp+1);    
    xc = xe(1:end-1) + dx/2;
    [xcm,ycm] = meshgrid(xc,xc);
    
    % Plot phi
    phi = reshape(amr_fp.data(2,:),mx_fp,my_fp)';
    pmin = min([pmin;phi(:)]);
    pmax = max([pmax;phi(:)]);
    hold on;
    if (min(phi(:)) < max(phi(:)))
        [~,h_phi(1)] = contour(xcm,ycm,phi,[0.5,0.5],'r','linewidth',1);
        hold on;
       %[~,h_phi(2)] = contour(xcm,ycm,phi,[tau/2, 1-tau/2],'k');
    end
end
    
    

fprintf('%10s %16.8e\n','pmin',pmin);
fprintf('%10s %16.8e\n','pmax',pmax);
fprintf('%10s %16.8e\n\n','pmax-pmin',qmax-qmin);
fprintf('%10s %16.8f\n','umin',qmin);
fprintf('%10s %16.8f\n','umax',qmax);
fprintf('%10s %16.8e\n','qmax-qmin',qmax-qmin);

showpatchborders;
setpatchborderprops('linewidth',1);
if (mx == 16)
    hidepatchborders(7);
elseif (mx == 32)
    hidepatchborders(8);
end

colormap(parula);

set(gca,'fontsize',16);

% create files for a movie;  Use 'convert2gif' to convert .jpg files
% to .gif files;  Use 'convert -delay 10 dendrite???.gif dendrite.gif'
% to create movie dendrite.gif.
NoQuery = 0;
MaxFrames = 64;
prt = true;
if (prt)
    if (length(amrdata) == 1)
        prefix = 'single';
    else
        prefix = 'amr10';
    end
    % hidepatchborders;
    title(sprintf('Thermal and phase fields (t = %.3f)',t),'fontsize',18);
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

shg

hold off;

clear afterframe;
