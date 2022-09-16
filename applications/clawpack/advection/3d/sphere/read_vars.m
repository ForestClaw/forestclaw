function parms = read_vars()

data = load('setprob.data');

parms.example = data(1);
parms.manifold = data(2);
parms.revs_per_sec = data(3);
parms.maxelev = data(4);

end