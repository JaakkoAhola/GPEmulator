# GPF

# Additions by Jaakko Ahola et al. 2020

This code is based on: [https://github.com/ots22/gpf/](https://github.com/ots22/gpf/)
We ran this code on Ubuntu 18.04.

1. First compile libraries. More information in section: Original Readme.

    - prerequisits
        1. install gfortran
        2. install cmake
        3. install build-essential

    - libraries are installed in a `$GPLIBRALY` of your choosing

    - LAPACK installation
        - `git clone https://github.com/Reference-LAPACK/lapack`
        - `cd lapack`
        - `cp make.inc.example make.inc`
        -  `ulimit -s unlimited` \# don't limit memory in compiling
        - `make`



    - NLOPT installation
        - `git clone https://github.com/stevengj/nlopt.git`
        - `cd nlopt`
        - `mkdir build`
        - `cd build`
        - `cmake ..`
        - `make`
        - `sudo make install`


    - makedepf90 installation
        - `sudo apt-get install makedepf90`

2. Installation of library files in this repository
    - `git clone https://github.com/JaakkoAhola/GPEmulator` \# use given or latest release
    - `make all`


3. Compile data reading library in the folder: [programs/lib_user](programs/lib_user)

4. Training the emulator is in folder: [programs/emulator_train](programs/emulator_train)
    * compile with command `make gp_train`
    * set input and output parameters in [Namelist](programs/emulator_train/train.nml)
    * run program with `./gp_train`
5. Predictor for emulator is in folder: [programs/emulator_predict](programs/emulator_predict)
    * compile with command: `make gp_predict`
    * set input and output parameters in [Namelist](programs/emulator_predict/predict.nml)

# Original readme:

GPF is a small Fortran library for Gaussian process regression.  It currently
implements value predictions with dense Gaussian processes, and projected-process
approximate Gaussian processes (described by [Rasmussen and Williams, 2006, chapter 8](http://www.gaussianprocess.org/gpml/chapters/RW8.pdf)).

## Installation (Linux)

GPF has the following external dependencies

* LAPACK
* NLopt (http://ab-initio.mit.edu/wiki/index.php/NLopt)
* makedepf90

GPF is written in standard Fortran 2008.  A makefile is provided for gfortran on linux.
It has been tested with gfortran 6.2, although earlier versions may also work.

Issuing
```sh
$ make all
```
will build the library (lib/libgpf.a) and tests.  The tests can be run with
```sh
$ make test
```

## Usage

### Examples

Some example programs using the library can be found in the directory `programs`.  See these
directories for further details.
* gp_train: train (and optionally optimize) a sparse GP from data files
* gp_predict: read in a GP file and produce outputs for points read from standard input

### Linking with the library

Ensure that the module files `gp.mod`, `gp_dense.mod` and `gp_sparse.mod` are in the
include path of your compiler.  Link with the static library `libgpf.a` by providing
the `-lgpf` flag (gfortran).

### Constructing a Gaussian process object from data

Include the module with either `use m_gp_dense` or `use m_gp_sparse`. These modules provide
the types `DenseGP` (full Gaussian process) and `SparseGP` (Gaussian process from the
projected process approximation).  Both are subtypes of the class `BaseGP`.

A `DenseGP` object can be constructed as follows:
```f90
type(DenseGP) my_gp
type(cov_sqexp) cf
type(noise_value_only) nm
my_gp = DenseGP(nu=[1.d-9], theta=[1.4_dp], x=x(1:N,:), obs_type=obs_type(1:N),
                t=t(1:N), CovFunction=cf, NoiseModel=nm)
```
where `nu`, `theta`, `x` and `t` are double precision arrays, and `obs_type` is an integer
array; `x` has rank two, the first dimension , and the second dimension defining the dimension of
the process.
* `nu`: the noise hyperparameters
* `theta`: the covariance hyperparameters
* `x`: the coordinates of the input data
* `obs_type`: the _types_ of the observations. If `obs_type(j)` is 0, then `t(j)` represents
an observation of the value of the underlying function.  If `obs_type(j)` is _i_ with _i > 0_,
then `t(j)` represents an observation of the partial derivative of the underlying function with
respect to the _i_ th component of _x_.
* `CovFunction`: the covariance function to use (of class `cov_fn`, see [below](#covariance-functions))
* `NoiseModel`: the noise model to use (of class `noise_model`, see [below](#noise-models))

A `SparseGP` object can be constructed similarly:
```f90
type(SparseGP) my_gp
my_gp = SparseGP(nsparse, nu, theta, x, obs_type, t, cf, nm)
```
* `nsparse`: an integer representing the number of privileged points to use in the projected
process. This should be strictly less than the number of elements in the training data set.
* The remaining parameters are the same as `DenseGP`.

### Reading a Gaussian process object from a file

```f90
DenseGP('in.gp')
```
or
```f90
SparseGP('in.gp')
```

### Writing a Gaussian process object to a file

```f90
class(BaseGP) my_gp
my_gp%write_out('out.gp')
```

### Maximum liklihood optimization of the hyperparameters

```f90
use gp_optim
class(BaseGP) my_gp
...
log_lik_optim(my_gp, lbounds, ubounds, optimize_max_iter, optimize_ftol)
```
* `my_gp`: an object of class BaseGP (either SparseGP or DenseGP)
* `lbounds` and `ubounds`: double precision arrays of length equal to the number of noise
hyperparemeters plus the number of covariance hyperparameters, representing lower and upper
bounds of the hyperparameters in the optimization.
* `optimize_max_iter`: the maximum number of iterations to perform in the optimization before
stopping.
* `optimize_ftol`: a double precision value stop the optimization when an optimization step
changes the log likelihood by less than this amount.

### Making predictions

A prediction of the underlying function at a coordinate `x` can be obtained as
```f90
y = my_gp%predict(x, obs_type)
```
* `x`: a double precision array with length of the dimension of the process
* `obs_type`: an integer representing the type of the observation. A value of `0` gives a
prediction of the value.  A value _i_ with _i > 1_, gives the predicted partial derivative of _y_
with respect to the _i_ th component of _x_.

### Covariance functions

Covariance functions extend the abstract type `cov_fn` (see [src/cov.f90]).

A square-exponential covariance function is defined as `cov_sqexp` ([src/cov_sqexp.f90]), and
others can be defined similarly.

### Noise models

Noise models extend the abstract type `noise_model`.

A noise model contains a function `noise`, taking the type of observation, the noise hyperparameters,
and returning the variance of the noise.

## References
Origal source codes: https://github.com/ots22/gpf/

C. Rasmussen and C. Williams. Gaussian Processes for Machine Learning. Adaptative
Computation and Machine Learning Series. MIT Press, 2006. ISBN 9780262182539.
