# ImpurityCIFbuilder

Often times we want to modify crystal structures with impurities (ex: doping). This script takes a crystalographic information file (cif) as input, uses it to create a new supercell, changes one of the atoms into an impurity, and then finds the new symmetry reduced cif. 

```
cif = "NiO_sym.cif"
supercell = [2,2,2]
impurity_elm = "Al"
impurity_coords = "center"
```

This is makes it easy to automate the generation of symmetry reduced cifs containing impurities which can then be transfered into any electronic structure code using something like cif2cell. This script uses the common libraries numpy and pymatgen.

## RUNNING

#### Impurity Element
The impurity element (`impurity_elm`) can be the same as one of the atoms already contained in the cif file. This is useful when you want a symmetry recuded version of the cif, but not something that is just a P1 spacegroup. This could be used when an electronic structure code only allows you to modify particular Wycoff positions, but you only want to change something about a single atom. 

For example

```
cif = "NiO_sym.cif"
supercell = [2,2,2]
impurity_elm = "Ni"
impurity_coords = "center"
```

gets you from a fully symmetrized cif with only 2 Wyckoff positions
```
Ni2+  Ni0  4  0.00000000  0.00000000  0.00000000  1 
O2-  O1  4  0.00000000  0.00000000  0.50000000  1
```

to a cif with 10 Wyckoff positions
```
Ni  Ni_imp  1  0.50000000  0.50000000  0.50000000  1
Ni2+  Ni0  12  0.00000000  0.25000000  0.25000000  1.0
Ni2+  Ni1  12  0.25000000  0.25000000  0.50000000  1.0
Ni2+  Ni2  3  0.00000000  0.00000000  0.50000000  1.0
Ni2+  Ni3  3  0.00000000  0.50000000  0.50000000  1.0
Ni2+  Ni4  1  0.00000000  0.00000000  0.00000000  1.0
O2-  O6  12  0.00000000  0.25000000  0.50000000  1.0
O2-  O7  8  0.25000000  0.25000000  0.25000000  1.0
O2-  O8  6  0.00000000  0.00000000  0.25000000  1.0
O2-  O9  6  0.25000000  0.50000000  0.50000000  1.0
```

#### Impurity Coordinates
The impurity coordinates (`impurity_coords`) can be either "center", "corner", or a 3 element list of abc fractional lattice coordinates. So "center" and [0.5,0.5,0.5] are equivalent definitions for the impurity coordinates.
