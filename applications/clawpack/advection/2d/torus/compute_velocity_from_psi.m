function [u,v,w] = compute_velocity_from_psi(psi_xi,psi_eta,tau1,tau2)
    
% Surface normal
nvec = cross(tau1,tau2);
nvec = nvec/norm(nvec,2);


% Compute grad_psi, where grad is a surface gradient
a11 = dot(tau1,tau1);
a22 = dot(tau2,tau2);
a12 = dot(tau1,tau2);

det = a11*a22 - a12*a12;
if det == 0
    error('velocity : Determinant is 0');
end

a11inv = a22/det;
a22inv = a11/det;
a12inv = -a12/det;
a21inv = -a12/det;
tau1inv = a11inv*tau1 + a12inv*tau2;
tau2inv = a21inv*tau1 + a22inv*tau2;

gradpsi = psi_xi*tau1inv + psi_eta*tau2inv;

% Compute velocity using  u = grad_psi x nvec
vel = cross(gradpsi,nvec);

u = vel(1);
v = vel(2);
w = vel(3);

end