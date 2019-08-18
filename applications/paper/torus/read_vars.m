function [example,ic,refine_pattern, alpha, beta,...
    revs_per_sec,cart_speed,theta_range,phi_range] = read_vars()

data = load('setprob.data');

example = data(1);
ic = data(2);
refine_pattern = data(3);
alpha = data(4);
beta = data(5);
revs_per_sec = data(6);
cart_speed = data(7);
theta_range = data(8:9);
phi_range = data(10:11);


end