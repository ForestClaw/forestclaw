function h = hsmooth(r,r0)

h1 = (tanh((r+r0)/0.01d0) + 1)/2.d0;
h2 = (tanh((r-r0)/0.01d0) + 1)/2.d0;

h = h1 - h2;

end

