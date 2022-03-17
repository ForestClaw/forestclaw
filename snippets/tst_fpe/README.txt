This example implements DYI exception handling.  It can be tested with

    $ tst_fpe <fpe_id>

where <fpe_id> is an FPE type to trap. 

Valid types are indicated by integer 1-5 : 

1 : INVALID, e.g. sqrt(-1), log(-1), acos(2) etc. 

2 : OVERFLOW, e.g. 2^1024, exp(1000)

3 : UNDERFLOW, e.g. denormalized numbers such as 2^(-1074), exp(-710)

4 : DIVBYZERO, e.g. 1/0

5 : Quiet NAN : This is useful for detecting use of uninitialized data

6 : Signaling NAN : This NAN will get trapped. 

Ideas for code taken from 

    https://developer.apple.com/forums/thread/689159

and

    https://stackoverflow.com/questions/69059981/how-to-trap-floating-point-exceptions-on-m1-macs





