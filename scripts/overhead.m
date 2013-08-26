function overhead(advance,exchange,regrid)

a = advance;
e = exchange;
r = regrid;

o = r + e;

fprintf('\n');
fprintf('%40s %6.1f%%\n','Over head as percentage of ADVANCE',100*o/a);
fprintf('%40s %6.1f%%\n','Over head as percentage of total',100*o/(a + o));
fprintf('\n');
fprintf(['Grid advances take %2.1f times as long as exchange and ',...
    'regridding\n'],a/o);

fprintf('\n');
