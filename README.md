# pc_multinest
A program combining Planck Likelihood Code (PLC) with the Cosmic Linear
Anisotropy Solving System (CLASS) to evaluate log likelihoods of different 
cosmological models.

MultiNest has been incorporated to scan parameter space and converge on the
most likely set of parameters.

## Contents of this README
- [Dependencies](#dependencies)
- [Branches](#branches)
- [Installation](#installation)
- [Running the Program](#running-the-program)
- [Troubleshooting](#troubleshooting)


## Dependencies
Following are dependencies for `pc_multinest`. Dependencies without a
hyperlink to their source are assumed to exist in some state on whichever
computing cluster that is being used. Those with links to their sources are
recommended to be compiled by hand.
- [CLASS](https://github.com/lesgourg/class_public)
  - OpenMP 
- [PLC](http://pla.esac.esa.int/pla/#cosmology)
  - [CFITSIO](http://heasarc.gsfc.nasa.gov/fitsio/)
  - LAPACK and BLAS, or Intel MKL
- [MultiNest](https://ccpforge.cse.rl.ac.uk/gf/project/multinest/)
  - OpenMPI


## Branches
Each branch of `pc_multinest` corresponds to a different configuration used on
a different computing cluster. Each branch also has a primary development focus:
- `master`: CoEPP cloud. Development primarily focussed on `pc_multinest`
- `local`: local machine. Development primarily focussed on `pc_speedtest`
- `phoenix`: Phoenix supercomputer. Development primarily focussed on
  `pc_multinest_mpi`


## Installation
> If you wish to make `pc_multinest` from scratch, you must first invent the
> Universe.
>
> &mdash; <cite>Carl Sagan</cite>

_DISCLAIMER: This installation guide was written primarily with Phoenix users_
_in mind. Do not be alarmed if whichever system you find yourself on is not_
_covered in this guide._

There are three dependencies for `pc_multinest`. We will run through how to
build and install each dependency properly such that `pc_multinest` will be able
to use them.

Before compiling any of the dependencies, it will first be necessary to clone
into `pc_multinest` using `git`. Run
```
$ git clone https://gitlab.coepp.org.au/harryp/pc_multinest.git
```
in your home directory to clone into this repository.

If you need to change to a particular branch, then run
```
$ git checkout <branch>
```

### CLASS
CLASS can be downloaded (or cloned) from the
[public GitHub repository](https://github.com/lesgourg/class_public). Assuming
`git` is already present on whichever computing cluster you are operating on,
to clone simply run
```
$ git clone https://github.com/hdp1213/class_public.git ~/class
```
which will create a new directory `class/` in the home directory containing the
CLASS code. This particular iteration of the code contains a branch with code
written to work with `pc_multinest`, specifically as it appears on Phoenix.

Before compilation, it is necessary to switch to the `phoenix` branch of CLASS:
```
$ git checkout phoenix
```

Then navigate to the `class` directory to run
```
$ make class
```
which will build CLASS. If this compiles without any errors, you're good to go!
Otherwise, you might want to examine the output and rectify it yourself.

### PLC
PLC can be downloaded in two parts from the
[ESA Planck Legacy Archive](http://pla.esac.esa.int/pla/#cosmology). Navigate to
the Likelihood tab and download both the _Planck_ Likelihood Code and the data
for the _Planck_ baseline results:
- `COM_Likelihood_Code-v2.0_R2.00.tar.bz2`
- `COM_Likelihood_Data-baseline_R2.00.tar.gz`

Extract both of these files using your uncompressor of choice to the home
directory. However, before building PLC it is first necessary to build CFITSIO.

#### CFITSIO
The CFITSIO source code can be downloaded from NASA's
[HEASARC Software page](http://heasarc.gsfc.nasa.gov/fitsio/). Extract the
`.tar` to your home directory. To build CFITSIO, run
```
$ ./configure
$ patch -p3 ~/cfitsio/Makefile < ~/pc_multinest/patches/cfitsio.patch
$ make shared
$ make install
```
where we have patched the CFITSIO `Makefile` after it was generated by the
`configure` script. This fixes an issue on Phoenix where the script initially
finds `gcc` instead of the Intel C compiler `icc`.

With CFITSIO installed, you can now run a similar build process on PLC.
Navigating back to the PLC directory, `plc-2.0/`, you should run
```
$ patch -p3 ~/plc-2.0/Makefile < ~/pc_multinest/patches/plc.patch
$ make
$ make install
```
to make PLC. Hopefully this should work without a hitch (depending on the
version of the compiler you are using). Otherwise, you will have to do a little
bit more to fix the problem. See the [Troubleshooting](#troubleshooting)
section at the end for more help.

Also be sure to peruse the PLC `Makefile` in case you want to specify or change
the location of any linear algebra packages you want to link against during
the build. The Intel packages should work for Phoenix after you have patched
the `Makefile`.

### MultiNest
Now you will need to make MultiNest. This can be downloaded from
[CCPForge](https://ccpforge.cse.rl.ac.uk/gf/project/multinest/), where you will
need to register for an account before downloading. Using an `.edu` address will
automatically validate the account so you don't have to wait for any manual
verification that you are currently associated in some regard to a
degree-granting institution.

After extraction of the files, you should run
```
$ patch -p3 ~/MultiNest_v3.10/Makefile < ~/pc_multinest/patches/multinest.patch
$ make
```

Here, the patch increases optimisation flags and changes the default `make`
target to the shared library object `libnest3.so` instead of the default
statically linked `libnest3.a`. It also turns on MPI support for Phoenix users.

### pc_multinest
At this point, you should finally be able to compile `pc_multinest`.
Congratulations! If this goes off without a hitch, you can consider yourself
very lucky and should immediately go out and buy a lottery ticket.

Before compilation, it is recommended to look at the `Makefile` to change the
root home directory for `pc_multinest`. You can also verify the other
directories correspond to where you have installed the other dependencies.

Depending on which branch you have checked out, it is wise to compile
whichever executable is the primary development focus. For the `phoenix` branch,
this is `pc_multinest_mpi`. So, by running
```
$ make pc_multinest_mpi
```
you should find you are the proud owner of a sparkly fresh copy of a
parallelised `pc_multinest`.

For the `master` branch, which is usually for use on the CoEPP cloud, you should
run
```
$ make pc_multinest
```

## Running the Program
`pc_multinest` takes some command line arguments used to specify where to write
MultiNest output files to. The full format is
```
./pc_multinest directory/to/output/root-
```

The first command line argument specifies the directory (`directory/to/output`)
and root (`root-`) to use for MultiNest output files. It defaults to
`output/pc_multinest-`. If you want to run multiple instances of `pc_multinest`
independently, you must change the root of the output files to something
different for each run, _e.g._ `output/pc_multinest_full_run00-`.

Giving no command line arguments to `pc_multinest` will just use default values
for each variable that can be changed. Right now, that is only the output
directory and root file name.

### Analysing Output
To analyse the MultiNest output from a successful run, `pc_multinest` makes use
of the [ROOT](http://root.cern.ch/) library and framework. An analysis script is
provided in the repo under the name `pc_plot.cc`. To run this, make sure the
`file_name` variable in the file points to the output file for the run you want
to analyse.

Output files in MultiNest are those with the form `<root>-.txt`. After checking
the file name in the script is correct, run the script by executing
```
$ root.exe -l pc_plot.cc
```

This command will boot up a ROOT instance, before executing the `pc_plot.cc`
script inside. If everything has gone well (and you have ROOT in the first
place), you should see the means and standard deviations of each of the
parameters MultiNest has scanned over.

#### `root.exe` not existing on Phoenix
Don't worry, it's most likely because you've forgotten to add the following to
your `~/.bashrc`:
```
module load ROOT/v6.07.06-intel-2015c-Python-2.7.10
```

This will load an Intel-compiled version of ROOT which will be Just That Tad
Bit Faster compared to the GCC version.

## Troubleshooting

### `<library>.so: cannot open shared object file`
If `pc_multinest` compiles but does not run because of an error of the form
```
./pc_multinest_mpi: error while loading shared libraries: <library>.so: cannot open shared object file: No such file or directory
```
it is because `<library>.so` is not on your `LD_LIBRARY_PATH` path environment
variable. In order to rectify this, you must find the location of the missing
`<library>.so` and add it to the path using
```bash
export LD_LIBRARY_PATH=/path/to/lib:$LD_LIBRARY_PATH
```

This can be put into your `.bashrc` so that it runs every time you open a new
`bash` session.

### PLC cannot find `clik_lensing.h` during compilation
This is a bit of a "bug" with PLC. You can remedy this by copying the
`clik_lensing.h` file from the `src/` directory into the `inc/` directory.
Compilation should then work.

### `forrtl: No such file or directory`
If you run `pc_multinest` only to find that MultiNest crashes the program with
an error of the form
```
forrtl: No such file or directory
forrtl: severe (29): file not found, unit 55, file </path/to/file>
```
it is because the containing directory for the file in question does not exist.
The problem can be rectified by making the containing directory and running the
program again.