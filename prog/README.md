# Extended Reference Interaction Site Model (XRISM) program: `rism0` #


`rism0` is an extended reference interaction site model (XRISM) program.
It computes various radial distribution functions $g(r)$
for simple molecular systems.



## Compiling ##

### Using Makefile ###

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

### Manual compiling ###

`rism0` can also be compiled manually:
```
gcc -O3 rism0.c -o rism0 -lfftw3 -lm
```
Or, if you do not have FFTW,
```
gcc -O3 -DNOFFTW rism0.c -o rism0_nofftw -lm
```


## Usage ##

### Basic usage ###

The main program is `rism0'.
The basic usage is
```
./rism0 your.cfg
```

Here, `your.cfg` is a configuration file for the model,
see `sample.cfg` for an example.
To customize a configuration file,
please use `rism.html` under the `web` directory.

You can also test a few stock models.
For example, to test model 14, type
```
./rism0 14
```
There are 16 stockmodels, and the last is chosen by default.

The default output is out.dat.
The output can be used in Gnuplot.



### Configure file ###

The configure file is a text file for the input.
It can be readily edited using any text editor,
such as `vim`, `emacs`, or `Notepad` under Windows.
The basic format is
```
key = value
```
where `key` specifies the variable name.

The best way to use a configuration file
is to modify an existing one,
e.g., `sample.cfg` in this directory.
The entries are quite self-explanatory.

Several points to remember

 * `key` or `value` names are case-insensitive.
 * Use `#` for comments.
```
# This is a comment line.
```
 * Arrays. `rho(1)` means the density of the first site.
 * Atom/site indices start from 1 instead of 0 (the first site is site 1).
 * The `data` directory contains many configuration files,
   many were used to reproduce existing literature.



### Override default parameters ###

Several commonly used parameters can be overridden on command line.
For example, if in the configuration file, the MDIIS solver is used,
and you want to change it to the LMV solver without editing the configuration file,
then type
```
./rism0 your.cfg -SLMV
```
The valid options are `Picard`, `LMV`, and `MDIIS`.

If you want to use the KH closure
```
./rism0 your.cfg -CKH
```
The valid options are `PY`, `HNC`, and `KH`.

If you want to change the charge of site 3 to 1.0 in your model
```
./rism0 your.cfg -q3,1.0
```
The site index starts from 1 instead of 0.

Similarly, if you want to change the density of site 2 to 0.03 in your model
```
./rism0 your.cfg -r2,0.03
```
The site index starts from 1 instead of 0.
If you want to add a bond between site 2 and site 3, with a distance 1.85
```
./rism0 your.cfg -d2,3,1.85
```
If you wish to annihilate an existing bond
```
./rism0 your.cfg -d2,3,inf
```

Type `./rism0 -h` for a full list of tunable parameters.



## Notes ##

### Reference implementations ###

Although `rism0` is program written from the scratch,
we have benefited from a similar program by Jesse Howard.

Pettitt Group internal reference:
  * /Bossman/Software/rism/
  * /Bossman/Software/3Drism/


### Water models ###

Currently, the data directory contains configuration files for the following water models
  * SPC
  * SPC/E
  * TIPS
  * TIP3P

External reference:
  * [Water models](https://en.wikipedia.org/wiki/Water_model)

Pettitt group internal link:
  * /Bossman/Software/3Drism/h2o_lib/


### On charged and/or polar systems ###

XRISM usually does not produce the correct $g(r)$
for charged and polar systems, or its logarithm,
the potential of mean force (PMF).
For this reason, there are several ways to the correct the PMF.
For various definitions of the PMF,
see the note under `doc/pmfs.dat`.


### On the handling of a small $g(r)$ ###

Under the HNC closure, $g(r) = exp[-v(r) + t(r)]$,
and the result is always nonnegative.

For a different PY/KH closure, the exact formula for $g(r)$ is $1 + c(r) + t(r)$,
which is not necessarily positive.
This can be inconvenient if we wish to compute the PMF from $-ln g(r)$.
Thus, if the output from first method $[1 + c(r) + t(r)]$ is negative or too small,
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
where the tenth column explicitly uses the $ln g(r)$.

In our test, the first method appears to produce more reasonable results than the second.
(Tested on Dec. 8, 2014.)


