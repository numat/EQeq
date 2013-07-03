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

Program expects two input files `ionization.dat` and `chargecenters.dat` as well as
a CIF file for an input structure. To run the program, pass the input CIF file path
as the first argument to the executable (i.e: `\EQeq_v1_00.exe MyDirectory/myfile.cif` )
Additional input parameters are optional. Please look at source code to see
what the other optional inputs are for (should be mostly self-explanatory).
