function parms = read_vars()

data = load('setprob.data');

parms.example = data(1);
parms.manifold = data(2);
parms.revs_per_second = data(3);
parms.longitude = data(4:5);
parms.latitude = data(6:7);

end