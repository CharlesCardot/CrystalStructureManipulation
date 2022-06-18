from pymatgen.io.cif import CifWriter
import numpy as np

elements = ["H","He","Li","Be","B","C","N","O","F","Ne"]

def find_unique_species(structure):
    """ Find unique species in a structure, returns a list of strings """
    
    species_in_struc = []
    for species in structure.species:
        san_species = ''.join([i for i in str(species) if i.isalpha()])
        if str(san_species) not in species_in_struc:
            species_in_struc.append(san_species)
    return species_in_struc


def make_impurity_cif(cif,structure,impurity_elm,impurity_coords,pos_tol=0.01):
    """ Builds and writes the impurity cif """
    
    # Used for impurities that are the same atomic symbol as the element being replaced
    ghost_impurity = ""
    unique_species = find_unique_species(structure)
    if impurity_elm in unique_species:
        for elm in elements:
            if elm not in unique_species:
                ghost_impurity = elm
                break

    imp_ind = -1
    if isinstance(impurity_coords,str):
        try:
            impurity_coords.isalpha()
        except ValueError:
            errmsg = "impurity_coords should only be defined as \"center\", \"corner\",\n"
            errmsg = errmsg+"or a list of [a,b,c] fractional lattice cordinates."
            print(errmsg)

        if impurity_coords.isalpha():
            if impurity_coords.lower() == "center":
                impurity_coords = [0.5,0.5,0.5]
            elif impurity_coords.lower() == "corner":
                impurity_coords = [0,0,0]
            else:
                raise ValueError("impurity_coords should only be defined as \"center\" or \"corner\" if written as a string\n")
    if isinstance(impurity_coords,list):
        impurity_coords = np.array(impurity_coords)
        for struc_key, atom_coords in enumerate(structure.frac_coords):
            max_coords = atom_coords*(1+pos_tol)
            min_coords = atom_coords*(1-pos_tol)
            if np.all(impurity_coords >= min_coords) and np.all(impurity_coords <= max_coords):
                imp_ind = struc_key
    else:
        raise TypeError("impurity_coords must be a string (\"center\", \"corner\") or a list of [a,b,c] fractional lattice cordinates.")

    if imp_ind == -1:
        errmsg = "impurity_coords " + str(impurity_coords) + " not found in structure.\n"
        errmsg = errmsg+"You must select the abc fractional lattice coords of one of the atoms\n"
        errmsg = errmsg+"in the new supercell to become the impurity."
        raise ValueError(errmsg)

    if ghost_impurity:
        structure[imp_ind] = ghost_impurity
    else:
        structure[imp_ind] = impurity_elm

    imp_cif = CifWriter(structure,symprec=0.1)
    filename = cif.replace(".cif","_"+impurity_elm+"_imp.cif")
    imp_cif.write_file(filename)
    imp_cif = open(filename,"r")
    imp_cif_lines = list(imp_cif.readlines())

    imp_line = ""
    imp_line_ind = -1

    # If a ghost impurity is used, replace it with the true impurity
    if ghost_impurity:
        for key, line in enumerate(imp_cif_lines):
            line_list = line.split()
            if line_list[0] == ghost_impurity:
                line_list[0] = impurity_elm
                line_list[1] = impurity_elm+"_imp"
                imp_line = "  " + "  ".join(line_list) + "\n"
                imp_line_ind = key
   
        # make the impurity the first atom in the file
        iter_var = imp_line_ind
        while imp_cif_lines[iter_var].strip()[0].isalpha():
            iter_var = iter_var - 1
        imp_cif_lines.insert(iter_var+1,imp_line)
        imp_cif_lines.pop(imp_line_ind+1)
        
    # Even if no ghost impurity is used, make the impurity the first atom in the file
    else:
        for key, line in enumerate(imp_cif_lines):
            line_list = line.split()
            if line_list[0] == impurity_elm:
                line_list[0] = impurity_elm+"_imp"
                imp_line = "  " + "  ".join(line_list) + "\n"
                imp_line_ind = key
        iter_var = imp_line_ind
        while imp_cif_lines[iter_var].strip()[0].isalpha():
            iter_var = iter_var - 1
        imp_cif_lines.insert(iter_var+1,imp_line)
        imp_cif_lines.pop(imp_line_ind+1)
   
    imp_cif.close()
    imp_cif = open(filename,"w")
    imp_cif.writelines(imp_cif_lines)
    imp_cif.close()
