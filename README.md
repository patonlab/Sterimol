Sterimol
=====

Python program for the calculation of [Sterimol](http://www.ccl.net/cca/software/SOURCES/FORTRAN/STERIMOL/) parameters: L, B1 and B5 for half-sandwich complexes and organic molecules. If used on half-sandwich complexes, it also generates [Tolman cone angles](https://en.wikipedia.org/wiki/Ligand_cone_angle) and metal to ring-centroid (unweighted) distances.

Developed by Dr Kelvin Jackson (Oxford) and [Dr Robert Paton](http://paton.chem.ox.ac.uk) (Oxford).



####Installation
1. Download the scripts from https://github.com/bobbypaton/Sterimol
2. Add the directory of the scripts to the PATH environmental variable (optional).  
3.	Run the script with Gaussian input or output files.

####Correct Usage

#####For half sandwich complexes

```python
sterimol.py file(s)
```
* This program will read Gaussian input or output files or half sandwich complexes.


#####For organic molecules

```python
sterimol.py (-a1 atom A) (-a2 atom B) (-radii radius-model) file(s)
```
* `-a1` and `-a2` specify atoms A and B atoms for the calculation - these fields are mandatory as they specify the axis along which Sterimol parameters are calculated.
* The `-radii` option specifies the radial model used. If left blank, the default model is the original CPK radii.


####Example 1:
Calculating Tolman cone angles, metal to ring-centroid distances, and Sterimol parameters for a half-sandwich complex from a Gaussian output file.

```python
./sterimol.py RhCpMe5Cl2PMe3.log

Sandwich Analysis
STERIMOL: using original CPK Van der Waals parameters

Structure                 Tolman_CA   MC_dist         L        B1        B5
RhCpMe5Cl2PMe3.com           173.97     1.833     4.016     3.902     4.304


```

The output shows the tolman cone angle (in degrees) and metal to centroid distance, L, B1 and B5 (all in Angstrom). Cone angles and Sterimol parameters are calculated using the original CPK atomic radii. 

####Example 2:
Calculating Sterimol parameters for a simple organic molecule (tert-butyl) from a Gaussian input.

```python
./sterimol.py -a1 2 -a2 1 tBu.com

STERIMOL: using original CPK Van der Waals parameters
Atoms 1 and 0 define the L-axis and direction [ 0.35667284  0.50439819 -0.8736515 ]

Atom       Xco/A     Yco/A     Zco/A    VdW/pm
##############################################
C         -0.893     1.083     0.000     150.0
H         -0.536     1.588    -0.874     100.0
C         -0.379     1.809     1.257     150.0
H          0.691     1.808     1.258     100.0
H         -0.734     2.819     1.256     100.0
H         -0.737     1.306     2.131     100.0
C         -2.433     1.084     0.000     150.0
H         -2.789     2.092     0.000     100.0
H         -2.789     0.579    -0.874     100.0
H         -2.789     0.579     0.873     100.0
C         -0.379    -0.368     0.000     150.0
H         -0.736    -0.873    -0.874     100.0
H          0.691    -0.368    -0.000     100.0
H         -0.736    -0.873     0.874     100.0

Structure                      L        B1        B5
tBu.com                     4.05      2.73      3.14
```

The output in this case returns the atom types, cartesian coordinates and atomic radii according to the CPK radial definitions. The sterimol parameters for the structure are given underneath; L, B1 and B5 are all given in Angstroms.

####Example 3:
Calculating parameters for a dimeric half-sandwich complex from a Gaussian output file.

```python
./sterimol.py Rh2Cl4IndStar2.log

Sandwich Analysis
STERIMOL: using original CPK Van der Waals parameters

Structure                 Tolman_CA   MC_dist         L        B1        B5
Rh17_dimer_wB97XD.log       191.283     1.763     6.184     3.381     5.607
Rh17_dimer_wB97XD.log       190.174     1.766     6.239     3.386     5.608


```

In this example two sets of parameters are produced - this occurs when the dimeric complex does not have a symmetry plane and thus measurements from each of the two metal centers yeilds different results. In the case of symmetric dimers, only a single set of parameters is generated (as they would be the same when measured from either metal center).


####Tips and Troubleshooting
* Errors will occur if this program is used on systems containing atoms for which there are no CPK defined radii.
* When running on organic molecules, the directionality of `-a1` and `-a2` is important - make sure the `-a1` to `-a2` vector is closer to the functional group being measured.
* It is possible to run on any number of files at once, for example using wildcards to specify all of the Gaussian files in a directory (*.out)
* The python file doesn’t need to be in the same folder as the Gaussian files. Just set the location of sterimol.py in the $PATH variable.




License: [CC-BY](https://creativecommons.org/licenses/by/3.0/)


