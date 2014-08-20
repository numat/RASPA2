RASPA2
======

A general purpose classical simulation package that can be used for the
simulation of molecules in gases, fluids, zeolites, aluminosilicates,
metal-organic frameworks, carbon nanotubes and external fields.

Installation
============

To download from Github and build from source, use:

```
git clone https://github.com/numat/RASPA2.git
cd RASPA2
mkdir -p m4
autoreconf --install
./configure --prefix=$PWD/simulations
make && make install
```

RASPA uses automake to manage code compilation, and the dependencies may not be
installed on your system by default. If this is the case, the install steps
will tell you on error. You can get the dependencies through your package
manager; the example below is for apt-get, but the command is similar for other
package managers (e.g. yum, brew).

```
sudo apt-get install automake autoconf libtool
```

Finally, in order for RASPA to find relevant files, you may need to configure
your environment variables. See below for more.

Usage
=====

There are currently two ways to use RASPA: through configuration files, or
through Python functions.

###Files

This approach uses RASPA script files (e.g. `simulation.input`), and a set of
directories to load structures, gases, and forcefield files. It outputs a set
of folders and files which contain data about the simulation run. This is the
original designed use of RASPA, and can be used in conjunction with shell
scripts for small sets of long-running simulations. To see the help, type:

```bash
/path/to/RASPA2/simulations/bin/./simulate -h
```

To write and configure simulation input files, read the documentation in the
`doc` folder.

###Python

The previous approach is useful for long-running jobs, but can become
cumbersome when you want to manage a large number of simulations. In this use
case, the Python wrapper provides a streamlined workflow.

RASPA's full functionality can be accessed in a Python script, which enables
simulation runs to be "glued" together with auxiliary workflow steps. For example,
to run a set of simulations across a logarithmic range of pressures, parse the
outputs for uptake data, and plot the results, use:

```python
import RASPA2

# Set up
gas = "CO2"
pressures = [1e4 * 10**(0.1 * i) for i in range(21)]

# Run
results = [RASPA2.run(my_structure, gas, temperature=298, pressure=pressure)
           for pressure in pressures]

# Parse
uptakes = [r["Number of molecules"][gas]
            ["Average loading absolute [cm^3 (STP)/cm^3 framework]"][0]
           for r in results]

# Plot
plot(pressures, uptakes)
```

When used in conjunction with supercomputer interfaces such as the
[IPython Notebook](http://ipython.org/notebook.html) and chemical data handling
libraries such as [open babel](http://openbabel.org/wiki/Main_Page), it is
possible to completely automate job distribution, cif formatting, charge
estimation, and more in a single script.

For more examples, check out the [NuMat workflow tutorial](https://github.com/numat/mofgen/wiki/Workflow)
or view some [example notebooks](https://github.com/numat/simulation-notebooks).

Environment
===========

Your environment is a set of variables that all programs can access, and can
be viewed by typing `env` into a terminal. Both RASPA and Python use specific
environment variables for operation, and you should ensure that the variables
mentioned are properly set up.

###RASPA2_DIR

By default, RASPA expects certain files to be in specific paths. Here's the logic:

1. See if that file is in your current directory.
2. See if that file is in a subdirectory of `RASPA2_DIR`.
3. See if that file is in a subdirectory of `RASPA_DIR`.
4. See if that file is in a subdirectory of `HOME`.

This is the case for everything without streaming, but it's also the case for
simulation gases and force fields with streaming. I like the solution:

```
cp -r simulations/share ~/
```

which moves the proper directory to `$HOME`, but you can also add this to your
`~/.bashrc`:

```
export RASPA2_DIR=/path/to/share
```

###PYTHONPATH

Python will look in certain directories when you try to import a library, e.g.
`import RASPA2`. The standard approach is to put your library in a directory
already in `PYTHONPATH`, such as:

```
mv RASPA2 ~/miniconda3/lib/python3.4/site-packages/RASPA2
```

(This example uses the [miniconda](http://conda.pydata.org/miniconda.html)
Python 3.4 distribution)

You can also set your own `PYTHONPATH` variable in your `~/.bashrc`:

```
export PYTHONPATH=/path/to/RASPA2
```

Be careful with this approach, though. If you don't keep your `~/.bashrc`
neatly organized, you will complicate your environment and make future installs
more painful.
