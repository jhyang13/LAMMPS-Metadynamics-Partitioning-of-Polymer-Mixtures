import os
import sys
import random
import math

natoms = int(sys.argv[1])
zlo = 0.858374
zhi = 120

def coord_str(atom_id, mol_id, atom_type, x, y, z):
    out_string = "%d %d %d %.6f %.6f %.6f\n"%(atom_id, mol_id, atom_type, x, y, z)
    return out_string

def bond_str(bond_id, bond_type, atom1_id, atom2_id):
        out_string = "%d %d %d %d\n"%(bond_id, bond_type, atom1_id, atom2_id)
        return out_string 

def angle_str(angle_id, angle_type, atom1_id, atom2_id, atom3_id):
				out_string = "%d %d %d %d %d\n"%(angle_id, angle_type, atom1_id, atom2_id, atom3_id)
				return out_string

def write_coords(natoms, zlo):
    out_string = ""
    for i in range(0,natoms,1):
        atom_id = i+1
        mol_id = 1 
        if i==0:
            atom_type = 1
        else:
            atom_type = 1
        x = 0.0
        y = 0.0
        z = zlo + i
        out_string += \
        coord_str(atom_id,mol_id,atom_type,x,y,z)
    return (out_string)

def write_bonds(natoms):
        out_string = ""
        for i in range(0,natoms-1,1):
                atom1_id = i+1
                atom2_id = i+2
                bond_id = i+1
                bond_type = 1
                out_string += \
                bond_str(bond_id,bond_type,atom1_id,atom2_id)
        return(out_string)

def write_angles(natoms):
				out_string = ""
				for i in range(0, natoms-2,1):
								atom1_id = i+1
								atom2_id = i+2
								atom3_id = i+3
								angle_id = i+1
								angle_type = 1
								out_string += angle_str(angle_id, angle_type, atom1_id, atom2_id, atom3_id)
				return(out_string)

out_string = "\n%d atoms\n"%(natoms)
out_string += "%d bonds\n"%(natoms-1)
out_string += "\n1 atom types\n"
out_string += "1 bond types\n"
out_string += "\n-15.00000 15.00000 xlo xhi\n"
out_string += "-15.00000 15.00000 ylo yhi\n"
out_string += "0.000000 %.6f zlo zhi\n"%(zhi)  
header_string = out_string                 

atom_return = write_coords(natoms, zlo)
out_string = "\nAtoms\n"
out_string +="  # atom_id, mol_id, atom_type, x, y, z \n"
atoms_string = header_string + out_string + atom_return

bond_return = write_bonds(natoms)
out_string = "\nBonds\n"
out_string +="  # bond_id, bond_type, atom1_id, atom2_id \n"
atoms_bonds_string = atoms_string + out_string + bond_return

angle_return = write_angles(natoms)
out_string = "\nAngles\n"
out_string += "  # angle_id angle_type atom1 atom2 atom3 \n"
atoms_bonds_angles_string = atoms_bonds_string + out_string + angle_return

f = open("chain.data","w+")
f.write(atoms_bonds_string)
f.close()
