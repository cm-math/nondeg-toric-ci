# nondeg-toric-ci
Magma programs used to perform computations in Theorem 1.3 from the
preprint "[The anticanonical complex
for non-degenerate toric complete intersections](https://arxiv.org/abs/2006.04723)" (arXiv:2006.04723v2)

## How to run the classification
Make sure you are in the directory with the file `classification.m`,
then start a Magma session and invoke `ClassifyTerminalFano3foldsFWPS`. 

```
> load "classification.m";  
Loading "classification.m"
> SD := ClassifyTerminalFano3foldsFWPS();
Start classification ...
```
While the compuations are running some output about
the number of verified candidates and the time taken (seconds) is printed.
```
x;u;nVerified;time
[ 1, 1, 1, 1, 1 ];[ 2 ];1;0.050
[ 1, 1, 1, 1, 1 ];[ 3 ];2;0.670
[ 1, 1, 1, 1, 1 ];[ 4 ];1;0.110
[ 1, 1, 1, 1, 2 ];[ 4 ];2;0.220
[ 1, 1, 1, 2, 2 ];[ 4 ];2;0.180
[ 1, 1, 1, 2, 2 ];[ 6 ];2;3.210
[ 1, 1, 1, 1, 3 ];[ 6 ];1;0.360
[ 1, 1, 1, 2, 3 ];[ 6 ];1;0.400
[ 1, 1, 2, 2, 3 ];[ 6 ];1;0.390
[ 1, 1, 2, 3, 3 ];[ 6 ];1;0.410
[ 1, 2, 2, 3, 3 ];[ 6 ];1;0.410
[ 1, 1, 1, 2, 4 ];[ 8 ];2;1.470
[ 1, 2, 3, 3, 4 ];[ 12 ];1;5.460
[ 1, 1, 3, 4, 4 ];[ 12 ];1;5.640
[ 1, 3, 3, 4, 4 ];[ 12 ];0;5.610
[ 1, 1, 2, 2, 5 ];[ 10 ];1;3.960
[ 1, 3, 3, 5, 5 ];[ 15 ];0;19.940
[ 1, 1, 2, 3, 6 ];[ 12 ];2;5.990
[ 1, 1, 1, 4, 6 ];[ 12 ];1;5.840
[ 1, 1, 3, 4, 6 ];[ 12 ];0;5.880
[ 1, 1, 2, 6, 9 ];[ 18 ];1;29.800
[ 1, 1, 4, 5, 10 ];[ 20 ];1;48.690
[ 1, 1, 3, 8, 12 ];[ 24 ];1;90.220
[ 1, 2, 3, 10, 15 ];[ 30 ];1;215.330
[ 1, 1, 6, 10, 15 ];[ 30 ];0;213.420
[ 1, 1, 6, 14, 21 ];[ 42 ];1;841.500
[ 1, 1, 1, 1, 1, 1 ];[ 2, 2 ];2;9.060
[ 1, 1, 1, 1, 1, 1 ];[ 3, 2 ];1;0.240
[ 1, 1, 1, 2, 2, 2 ];[ 4, 4 ];4;20.270
[ 1, 2, 2, 2, 3, 3 ];[ 6, 6 ];1;7.170
[ 1, 1, 2, 3, 3, 3 ];[ 6, 6 ];1;6.940
[ 1, 2, 2, 3, 3, 3 ];[ 6, 6 ];1;7.190
[ 2, 2, 2, 3, 3, 3 ];[ 6, 6 ];0;6.110
[ 2, 2, 2, 5, 5, 5 ];[ 10, 10 ];0;81.080
[ 1, 1, 1, 1, 1, 1, 1 ];[ 2, 2, 2 ];4;4446.860
```
Finally, the function returns a list of specifying data,
i.e., triples (Q, mu, T) where
* Q = [w_1, .., w_r] is the (generator) degree matrix,
* mu = [mu_1, ..., mu_s] is the relation degree matrix,
* T = [t_1, ..., t_q] contains the orders of the finite cyclic factors of K = Cl(Z)

```
> #SD;
42
```
We examine some arbitrary entries.
```
> SD[1];
[*
    [
        [ 1, 1, 1, 1, 1 ]
    ],
    [
        [ 2 ]
    ],
    []
*]

> SD[3];
[*
    [
        [ 1, 1, 1, 1, 1 ],
        [ 0, 0, 1, 1, 2 ]
    ],
    [
        [ 3, 0 ]
    ],
    [ 3 ]
*]

> SD[42];
[*
    [
        [ 1, 1, 1, 1, 1, 1, 1 ],
        [ 0, 0, 0, 0, 1, 1, 1 ],
        [ 0, 0, 1, 1, 0, 0, 1 ],
        [ 0, 1, 0, 1, 0, 1, 0 ]
    ],
    [
        [ 2, 0, 0, 0 ],
        [ 2, 0, 0, 0 ],
        [ 2, 0, 0, 0 ]
    ],
    [ 2, 2, 2 ]
*]
```
