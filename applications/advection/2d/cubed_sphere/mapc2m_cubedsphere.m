function [xp,yp,zp] = mapc2m_cubed_sphere(xc,yc)

% static inline void
% fclaw2d_map_c2m_csphere_help (double R, double xi, double eta,
%                              double *x, double *y, double *z)
% {
%     const double tan_xi = tan (.5 * M_PI * (xi - .5));
%     const double tan_eta = tan (.5 * M_PI * (eta - .5));
% 
%     FCLAW_ASSERT (0. <= xi && xi <= 1.);
%     FCLAW_ASSERT (0. <= eta && eta <= 1.);
% 
%     *z = R / sqrt (SC_SQR (tan_xi) + SC_SQR (tan_eta) + 1.);
%     *x = *z * tan_xi;
%     *y = *z * tan_eta;
% }


blockno = getblocknumber();

switch blockno
    case 0
        [yp,xp,zp] = csphere_basic(xc,yc);
        zp = -zp;
    case 1
        [zp,xp,yp] = csphere_basic(xc,yc);
    case 2
        [zp,yp,xp] = csphere_basic(xc,yc);
        xp = -xp;
    case 3
        [xp,yp,zp] = csphere_basic(xc,yc);
    case 4
        [xp,zp,yp] = csphere_basic(xc,yc);
        yp = -yp;
    case 5
        [yp,zp,xp] = csphere_basic(xc,yc);
    otherwise
end
                

end

function [xp,yp,zp] = csphere_basic(xc,yc)

tan_xi = tan(0.5*pi*(xc-0.5));
tan_eta = tan(0.5*pi*(yc-0.5));
zp = 1./sqrt(tan_xi.^2 + tan_eta.^2 + 1);
xp = zp.*tan_xi;
yp = zp.*tan_eta;

end
