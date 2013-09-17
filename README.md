RASPA2
======

A general purpose classical simulation package that can be used for the
simulation of molecules in gases, fluids, zeolites, aluminosilicates,
metal-organic frameworks, carbon nanotubes and external fields. For more, see
the `README`.

Installation
============

This can be done through the standard automake procedure, with some small
edits.

```bash
cd /path/to/RASPA2
cp src/Makefile.am.object src/Makefile.am
mkdir -p m4
autoreconf --install
./configure --prefix=$PWD/simulations
make && make check && make install
```

You then should set a `RASPA_DIR` environment variable. In the case of the
above install directory...

```bash
export RASPA2_DIR=/path/to/RASPA/simulations
```

You can now run RASPA2 simulations through the installed binary. Check the help.

```bash
/path/to/RASPA2/simulations/bin/./simulate -h
```

Streaming + Shell Scripting
===========================

In this version of RASPA2, you can run everything as a stream, allowing
inclusion into a bash script without dealing with files and folders. The
above command could be streamed with the `-s` flag, replacing input files
with data.

```bash
./simulate -s -i "`cat test.input`" -c "`cat mof.cif`" > output.txt
```

While this looks more complicated, you can write a bash script using variables
without requiring filesystem commands.

```bash
script=`cat test.input`
mof=`cat mof.cif`
output=`./simulate -s -i "$script" -c "$mof"`
```

Python Bindings
===============

For use as one step in a high-throughput process, RASPA2 can be compiled as a
shared object and interfaced with other programs through Python "glue". To do
this,

```bash
cd /path/to/RASPA2
cp src/Makefile.am.shared src/Makefile.am
mkdir -p m4
autoreconf --install
./configure --prefix=$PWD/simulations
make && make check && make install
sudo cp simulations/lib/libraspa2.so /usr/lib
```

(if you don't have sudo access, `mkdir ~/lib; cp libraspa2.so ~/lib`, then
change the path at the top of `raspa2.py`)

This creates `libraspa2.so`, which is accessible to external programs. To use
the Python program via the command line, try

```
python raspa2.py *input_file* *structure_file*
```

This API allows for the streaming of file output, instead letting the user pipe
their results into whatever they see fit. This enables databases, sockets, and
a range of other possibilities.

The Python bindings can also be used as part of a larger "glue" script. To do
this, use

```python
import RASPA2
results = RASPA2.run(my_structure, "CO2", temperature=298.0, pressure=1e5,
                    helium_void_fraction=0.45)
```

This allows streaming without any required input or output files, which makes
it much easier to use with modern supercomputers (ie. [PiCloud](https://www.picloud.com))
and modern interfaces (ie. [IPython Notebook](http://ipython.org/notebook.html)).

RASPA2_DIR
=========

######Note: For this version, I've renamed `RASPA_DIR` to `RASPA2_DIR`. This allows for RASPA and RASPA2 to peacefully coexist. When RASPA is completely deprecated, I may fix the nomenclature.

By default, RASPA2 expects certain files to be in specific paths. Here's the logic:

1. See if that file is in your current directory `$PWD`.
2. See if that file is in a subdirectory of `$RASPA2_DIR`.
3. See if that file is in a subdirectory of `$HOME`.

This is the case for everything without streaming, but it's also the case for
simulation gases and force fields with streaming. I like the solution:

```
cp -r simulations/share ~/
```

but you can also add this to your `.bashrc`:

```
export RASPA2_DIR=/path/to/share
```

or pass an argument to the compiled code.

```
./simulate -d /path/to/share
```

(The latter functionality is also in the Python bindings, but I'm hiding it for
the sake of simplicity. If you *really* want to use it this way, just tell me.)
