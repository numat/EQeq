EQeq
====

Charge equilibration method for crystal structures

The source code in this program demonstrates the charge equilibration method described
in the accompanying paper. The purpose of the source code provided is to be
minimalistic and do "just the job" described. In practice, you may wish to add various
features to the source code to fit the particular needs of your project.

#####Major highlights of program:

 * Obtains charges for atoms in periodic systems without iteration
 * Can use non-neutral charge centers for more accurate point charges
 * Designed for speed (but without significant code optimizations)

#####Features not implemented but that you may want to consider adding:

 * Spherical cut-offs (for both real-space and reciprocal-space sums)
 * An iterative loop that guesses the appropriate charge center (so the user does not have to guess)
 * Ewald parameter auto-optimization
 * Various code optimizations

#####Running the program:

Program expects two input files `ionization.dat` and `chargecenters.dat`. Please
look at source code to see what the other optional inputs are for (should be
mostly self-explanatory). Compile with something like:

```
g++ main.cpp -O3 -o eqeq
```

and run with

```
./eqeq my_file.cif
```

#####Python bindings (BETA)

To facilitate automation and scaling, this version of EQeq can be operated via
a Python function. To enable, you must build EQeq as a shared library:

```
g++ -c -fPIC main.cpp -o eqeq.o

# Linux
g++ -shared -Wl,-soname,libeqeq.so -o libeqeq.so eqeq.o
# Mac
g++ -shared -Wl,-install_name,libeqeq.so -o libeqeq.so eqeq.o
```

Then, you must put the `EQeq` directory on your PYTHONPATH (remember to put this
in your ~/.bashrc to make it permanent).

```
export PYTHONPATH:/path/to/EQeq:$PYTHONPATH
```

Then, you can call EQeq from Python with

```python
import EQeq

EQeq.run("IRMOF-1.cif")
```

The input takes both filenames and actual data. The latter is useful when
connecting EQeq with other programs (which I plan to do eventually).
