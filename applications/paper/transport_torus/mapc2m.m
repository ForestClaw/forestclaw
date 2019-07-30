function [xp,yp,zp] = mapc2m(xc,yc)

% d = load('setprob.data');
% alpha = d(2);
% beta = d(3);

alpha = 0.4;
beta = 0.5;

[xp,yp,zp] = mapc2m_torus(xc,yc,alpha,beta);
   
end
