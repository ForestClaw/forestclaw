function uo = underover()

% UNDEROVER returns structure for assigning colors to under/overshoots
%
% UO = UNDEROVER().  To have special highlighting for under/over shoot
% values, the user should first set the ShowUnderOverShoots flag in
% SETPLOT2 or the SETPLOT3 files.
% Then, make a local copy of this file and supply values for the the
% following fields in the UO structure :
%
%          color_under     the rgb triple used to color undershoots
%          color_over      the rgb triple used to color overshoots
%          value_lower     lower value of exact solution
%          value_upper     upper value of exact solution
%          tol             tolerance for under/over shoots
%          colormap        colormap for values in [qlow-tol,qhi+tol];
%
% Values are considered in the overshoot/undershoot region if they satisfy
%
%                 q < value_lower - tol      (undershoot)
% or
%                 q > value_upper + tol      (overshoot)
%
% Values not in the under/overshoot region are linearly scaled into the
% colormap.  Values that lie in [qlow-tol,qlow] are clamped to qlow, and
% values in [qhi, qhi + tol] are clamped to qhi.
%
% Example :
%
%            function uo = underover()
%
%            colormap('default');
%            cm = colormap;
%
%            uo = struct('color_under',[0 1 1],...  % cyan
%                        'color_over',[1 0 1],...   % magenta
%                        'value_lower',0, ...
%            	         'value_upper',1,...
%            	         'tol',0,...
%            	         'colormap',cm);
%
% This function is called from SETCOLORS if the ShowOverUnderShoots flag is
% set.
%
% See also SETCOLORS, UNDEROVER_COLORBAR.
%

cm = colormap;

uo = struct('color_under',[0 1 1],...
            'color_over',[1 0 1],...
            'color_nan',[1 1 1], ...
            'value_lower',-1, ...
    	    'value_upper',1,...
            'tol',1e-5,...
            'colormap',cm);
