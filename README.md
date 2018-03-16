***********************
## DINAMICA
***********************

Dinamica is intended to provide a general interface to the analysis of the
Dynamical Systems. However, the specific set of aims ranges from the determining of
the dynamical regime of the Deterministic (Coupled Symmetric) Dynamical systems
written in the ODE form to the comparative analysis of the deterministic and
stochastic (both discrete and continuous) representations of the same system. The
aims are achieved by calculating various errors and comparing them with
user-defined levels. Thus, the user ought to have some experience in scientific
terms, although the algorithms are straightforward.

Dinamica treats both the stationary and oscillatory dynamical regimes.
Additionally, there are algorithmic ways to calculate various oscillatory regimes:
simple self-oscillations, coupled synchronized and non-synchronized, homogeneous
and non-homogenous regimes, to name a few.

Additionally, there are many techniques used assisting the user in the analysis of
dynamical systems. For example, period calculation with several methods, lyapunov
exponents, continuation over parameters of the system, different statistical tools,
period distributions etc.

Dinamica uses the powerful scientific library GSL (GNU Scientific Library) for the
routine operations like e.g. numerical integration, random number generation,
finding max and min etc. So one needs the GSL library pre-installed in order to
use Dinamica.

**********************
### Pre-requisites
**********************

1. Unix/Linux OS.
2. gcc compatible C compiler, needed for Dinamica functioning (not only
compilation).
3. GNU Scientific Library (GSL) installed.
4. Gnuplot plotting utility installed (optional, but advisable).

Quick installation:
1. Dowload the archive (usually .zip or .tar.gz) and uncompress it.
1a. NOTE: this repository contains the configuration files generated with ver. 1-14 of the GNU autotools.
If your version differs you need to run `aclocal` and `automake` first to regenerate the
configuration files for the updated autotools.
2. Configure the systemm by typing `./configure`. This will check for all the
requirements and complain if any of those is not found. You may consider `CPPFLAGS`
and `LDFLAGS` variables, as well as `--prefix` option to `./configure`, before
configuring the system (see below).
3. Type `make` to compile the `libdin.a` and the dinamica itself. The two must
appear under the `src/` directory in the root.
4. Type `make install` to install dinamica executable and `libdin.a` library to the
usual destinations (`/usr/local/bin` and `/usr/local/lib`, respectively). This might
require the root password. You may uninstall the program later by typing `make
uninstall` to remove those two files from the system. After the installation one
might want to remove all the files extracted from the archive.


********************************************************
### Important to know before configuring the package
********************************************************
Dinamica processes the input from the user (equations, parameters, variables,
constants etc.) in the form of the *.ode* file. The result of the processing is the
output *.c* file with C language definitions and functions for the user system and
the binary configuration *.bcf* file. This output *.c* file is then compiled with the
dinamica library (libdin.\*) generating the final executable *.din*. This executable
is then invoked to read the *.bcf* file and, finally, the program fires up.

It is important to understand that the compilation and linking of the libraries
take place during the functioning of Dinamica. Thus, the compiler and the right
path for the required libraries are needed to be properly set when Dinamica is
first compiled. One should take care of this before the configuration starts.

The way Dinamica can obtain the full set of the paths is to specify `CPPFLAGS`,
`LDFLAGS` and `--prefix` option, together or one by one when needed. These three
variables are passed to the Dinamica compilation command, so they are crucial. The
current directory where Dinamica is invoked is always checked for the dinamica
library (libdin.\*).


#### CPPFLAGS
This variable is important for the preprocessor, a program checking the included,
so called header, files. During the configuring the `./configure` script checks for
several *.h* files from the GSL and the standard C libraries whether they are
available. Sometimes it fails to find them in the standard locations. In this case,
the `CPPFLAGS` is needed. The usage is simple: if one knows that the GSL headers are
located, for example, in `/usr/local/include/`, e.g. the full path to `gsl_odeiv2.h`
file is `/usr/local/include/gsl/gsl_odeiv2.h` (similarly for all other gsl *\*.h*
files), then one could type:

```CPPFLAGS=-I/usr/local/include ./configure```

in the shell prompt. This sets the environment variable for the `./configure` script.
Next, in Dinamica this option will be set, i.e.  will be used as a path to the header
files to find. Note, that the path to the GSL header files is always *gsl/\*.h*, so omit
the `gsl/` directory when specifying the CPPFLAGS.

#### LDFLAGS
This variable shows the path to the GSL library (and possibly the standard C
libraries). If the gsl library files are located in `/usr/local/lib` then combining
with `CPPFLAGS` one could type

```CPPFLAGS=-I/usr/local/include LDFLAGS=-L/usr/local/lib ./configure```

which should do the trick. Additionally, Dinamica will compile the output *.c* file
using these options or the one specified by `--prefix` (see below) to find the
`libdin.*` library.


#### --prefix
`./configure` script accepts the `--prefix` option, which specifies the installation
directory for `make install` command. Namely,

```./configure --prefix=/usr/local```

would install all the files produced by the package to the corresponding
subdirectories of `/usr/local`. For example, binary files, i.e dinamica, would go
into `/usr/local/bin` and libraries, i.e. `libdin.a`, --- into `/usr/local/lib`. This
variable is set after the `-L` flag to the compiler in Dinamica.


*****************
### ODE FILES
*****************
The main input from the user of Dinamica is represented by the ode-files. The
ode-file syntax in many ways resembles that of the XPPAUT *.ode* files (Bard
Ermentrout's program for exploring the dynamical systems, see
<http://www.math.pitt.edu/~bard/xpp/xpp.html>). However, there are many Dinamica
specific options, like, for example, `#system` directive which instructs Dinamica how
many physical coupled systems are there. See more on how to write ode-files in the
manual. To launch the program type:

```dinamica <your_file>.ode```

this should process the file, compile and fire up Dinamica. There are several
command line options and many other things one can do with Dinamica. For those,
consult the manual

**************************************************************
Bug reports, any questions and suggestions can be sent to:
elias.potapov@gmail.com
**************************************************************
