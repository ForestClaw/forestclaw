function [example, beta_choice, alpha, x0, y0, a,b,...
    eps_disk,m_polar,x0_polar,y0_polar,r0_polar,r1_polar,...
    n_polar, bc] = read_vars()

data = load('setprob.data');

example = data(1);
beta_choice = data(2);
alpha = data(3);
x0 = data(4);
y0 = data(5);
a = data(6);
b = data(7);

% Polar data
eps_disk = data(8);
m_polar = data(9);
m = m_polar;
x0_polar = data(10+0*m:10+1*m-1); 
y0_polar = data(10+1*m:10+2*m-1);
r0_polar = data(10+2*m:10+3*m-1);
r1_polar = data(10+3*m:10+4*m-1);
n_polar  = data(10+4*m:10+5*m-1);

bc = data(10+5*m:end);


end