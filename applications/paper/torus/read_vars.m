function [example,ic,refine_pattern, alpha, beta,init_radius, ...
    revs_per_sec,cart_speed,theta_range,phi_range] = read_vars()

data = load('setprob.data');

example = data(1);
ic = data(2);
refine_pattern = data(3);
alpha = data(4);
beta = data(5);
init_radius = data(6);
revs_per_sec = data(7);
cart_speed = data(8);
theta_range = data(9:10);
phi_range = data(11:12);


end