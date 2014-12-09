# Extended Reference Interaction Site Model (XRISM) program: `rism0` #


`rism0` is an extended reference interaction site model (XRISM) program.
It computes various radial distribution functions g(r)
for simple molecular systems.



## Compiling ##

Open a terminal and enter
```
make
```
This produces an executable `rism0`.

The default `make` assumes that you have FFTW3 (`-lfftw3`) installed.
If this is not so, try
```
make nofftw
```
This produces `rism0_nofftw`.



## Basic usage ##

The main program is `rism0'.
The basic usage is
```
./rism0 your.cfg
```

Here, `your.cfg` is a configuration file for the model,
see `sample.cfg` for an example.
To customize a configuration file,
please use `gencfg.html` under the `web` directory.

You can also test a few stock models.
For example, to test model 14, type
```
./rism0 14
```
There are 16 stockmodels, and the last is chosen by default.

The default output is out.dat.
The output can be used in Gnuplot.



## Notes ##

### On charged and/or polar systems ###

XRISM usually does not produce the correct g(r)
for charged and polar systems, or its logarithm,
the potential of mean force (PMF).
For this reason, there are several ways to the correct the PMF.
For various definitions of the PMF,
see the note under `doc/pmfs.dat`.


### On the handling of a small g(r) ###

Under the HNC closure, g(r) = exp[-v(r) + t(r)],
and the result is always nonnegative.

For a different PY/KH closure, the exact formula for g(r) is 1 + c(r) + t(r),
which is not necessarily positive.
This can be inconvenient if we wish to compute the PMF from -ln g(r).
Thus, if the output from first method [1 + c(r) + t(r)] is negative or too small,
we automatically switch to the HNC formula to guarantee a nonnegative output.

This expedient practice appears to work for practical purposes.
A test case (using stock model 16) is the following.

```
make && ./rism0 16 -C KH -o out.dat
gcc -DGR_TINY=0 -O3 rism0.c -lfftw3 -lm -o rism0b && ./rism0b 16 -C KH -o out0.dat
```

Here, on the first line, we compile the program `rism0` with the default setting,
then run it on model 16.
The option `-C KH` asks the program to use the KH closure.
On the second line, we compile a different `rism0b`.
The compilation flag `-DGR_TINY=0` manually sets the threshold to 0.0
instead of the default value `100 * model->tol`.

For a comparison, enter Gnuplot, and type
```
plot [:10][-20:20] "out.dat" u 1:10 w l, "out0.dat" u 1:10 w l
```
where the tenth column explicitly uses the ln g(r).

In our test, the first method appears to produce more reasonable results than the second.
(Tested on Dec. 8, 2014.)
