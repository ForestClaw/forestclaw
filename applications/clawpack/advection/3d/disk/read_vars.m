function parms = read_vars()

data = load('setprob.data');

parms.example = data(1);
parms.manifold = data(2);
parms.alpha = data(3);

end