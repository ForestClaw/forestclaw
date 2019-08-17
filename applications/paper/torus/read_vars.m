function [example,ic,refine_pattern, alpha, beta,revs_per_sec] = read_vars()

data = load('setprob.data');

example = data(1);
ic = data(2);
refine_pattern = data(3);
alpha = data(4);
beta = data(5);
revs_per_sec = data(6);


end