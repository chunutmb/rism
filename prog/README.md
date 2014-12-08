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



## Notes on charged and/or polar systems ##

XRISM usually does not produce the correct g(r)
for charged and polar systems, or its logarithm,
the potential of mean force (PMF).
For this reason, there are several ways to the correct the PMF.
For various definitions of the PMF,
see the note under `doc/pmfs.dat`.

