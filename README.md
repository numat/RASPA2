RASPA2
======

A general purpose classical simulation package that can be used for the
simulation of molecules in gases, fluids, zeolites, aluminosilicates,
metal-organic frameworks, carbon nanotubes and external fields.

Installation
============

RASPA uses automake to manage code compilation. There are currently two separate
ways to run RASPA, which each have slightly different build procedures.

###Files

This approach uses RASPA script files (e.g. `simulation.input`), and a set of
directories to load structures, gases, and forcefield files. It outputs a set
of folders and files which contain data about the simulation run. This is the
original designed use of RASPA, and can be compiled like so:

```
cd /path/to/RASPA2
cp src/Makefile.am.object src/Makefile.am
mkdir -p m4
autoreconf --install
./configure --prefix=$PWD/simulations
make && make check && make install
```

You will also want to handle the `RASPA2_DIR` environment variable. See below for more.

Once the variable is set, you can run RASPA2 simulations through the installed binary. Check the help.

```bash
/path/to/RASPA2/simulations/bin/./simulate -h
```

###Python

Instead of loading and writing files, this approach uses Python functions to
run simulations. This is especially useful in high-throughput screening,
where the goal is to compile the results of many separate simulations into a
final output. To build:

```
cd /path/to/RASPA2
cp src/Makefile.am.shared src/Makefile.am
mkdir -p m4
autoreconf --install
./configure --prefix=$PWD/simulations
make && make check && make install
```

You will want to handle the `RASPA2_DIR` environment variable, as well as
manage your `PYTHONPATH`. See below for more.

This can still be run from the command line (see `python raspa2.py -h` for help),
but is intended for use as part of a Python script in e.g. an
[IPython Notebook](http://ipython.org/notebook.html). For example, to run a set
of simulations across a logarithmic range of pressures, parse the outputs for
uptake data, and plot the results, use:

```python
import RASPA2

# Set up
gas = "CO2"
temperature = 298
pressures = [1e4 * 10**(0.1 * i) for i in range(21)]
results = []

# Run
for pressure in pressures:
    results.append(RASPA2.run(my_structure, gas, temperature=temperature,
                              pressure=pressure))

# Parse
uptakes = [r["Number of molecules"][gas]
            ["Average loading absolute [cm^3 (STP)/cm^3 framework]"][0]
           for r in results]

# Plot
plot(pressures, uptakes)
```

For more examples, check out the [NuMat workflow tutorial](https://github.com/numat/mofgen/wiki/Workflow)
or view some [example notebooks](https://github.com/numat/simulation-notebooks).

Environment
===========

Your environment is a set of variables that all programs can access. You
can view it by typing `env` into a terminal. Both RASPA and Python use specific
environment variables for operation, and you should ensure that the variables
mentioned are properly set up.

###RASPA2_DIR

By default, RASPA expects certain files to be in specific paths. Here's the logic:

1. See if that file is in your current directory `$PWD`.
2. See if that file is in a subdirectory of `$RASPA2_DIR`.
3. See if that file is in a subdirectory of `$RASPA_DIR`.
4. See if that file is in a subdirectory of `$HOME`.

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
