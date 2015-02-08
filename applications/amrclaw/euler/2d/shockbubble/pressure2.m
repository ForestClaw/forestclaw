function pressure = pressure2(data)
% compute the pressure from the data in 2 dimensions

gamma = 1.4;
rho = data(:,1);
u = data(:,2)./rho;
v = data(:,3)./rho;
energy = data(:,4);
pressure = (gamma-1) * (energy - 0.5*rho.*(u.*u + v.*v));
