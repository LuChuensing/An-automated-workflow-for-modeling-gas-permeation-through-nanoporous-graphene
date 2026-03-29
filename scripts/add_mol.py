import numpy as np
from ase.io import read, write
from ase import Atom
import argparse
import sys

def find_hole_center_by_edge_H(graphene, lx, ly):
    h_positions = np.array([atom.position for atom in graphene if atom.symbol == 'H'])
    if len(h_positions) == 0: 
        return lx/2, ly/2
    
    # Calculate center via periodic mean
    theta_x = (h_positions[:, 0] / lx) * 2 * np.pi
    theta_y = (h_positions[:, 1] / ly) * 2 * np.pi
    center_x = (np.arctan2(np.mean(np.sin(theta_x)), np.mean(np.cos(theta_x))) / (2 * np.pi) * lx) % lx
    center_y = (np.arctan2(np.mean(np.sin(theta_y)), np.mean(np.cos(theta_y))) / (2 * np.pi) * ly) % ly
    return center_x, center_y

def main():
    parser = argparse.ArgumentParser(description="Smart Gas Insertion Tool")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-g", "--gas", required=True, choices=['H2', 'N2', 'CH4', 'Ar', 'CO2'])
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    system = read(args.input)
    lx, ly, _ = system.cell.lengths()
    cx, cy = find_hole_center_by_edge_H(system, lx, ly)
    cz = 65.0 

    n_graphene_atoms = len(system)

    # -------------------------------------------------------------------
    # SMART FIX: Proxy Elements
    # We use 'F' and 'Cl' to prevent ASE from merging gas atoms with 
    # Graphene atoms (like the 'C' in CO2 clashing with Graphene 'C').
    # LAMMPS only reads the integer Types (1, 2, 3, 4) anyway.
    # -------------------------------------------------------------------
    sym_outer = 'F'
    sym_inner = 'Cl'

    # Dictionary: [bond_length, is_3_site]
    gas_params = {
        'H2':  [0.74, True], 
        'N2':  [1.10, True], 
        'CO2': [1.16, True], 
        'CH4': [0.0,  False],  
        'Ar':  [0.0,  False]   
    }

    bond_len, is_3_site = gas_params[args.gas]

    if is_3_site:
        if args.gas == 'CO2':
            # Linear Outer-Inner-Outer (O-C-O proxy)
            system.append(Atom(sym_outer, position=(cx, cy, cz - bond_len)))
            system.append(Atom(sym_outer, position=(cx, cy, cz + bond_len)))
            system.append(Atom(sym_inner, position=(cx, cy, cz)))
        else:
            # Linear Outer-Inner-Outer (X-COM-X proxy)
            r = bond_len / 2.0
            system.append(Atom(sym_outer, position=(cx, cy, cz - r)))
            system.append(Atom(sym_outer, position=(cx, cy, cz + r)))
            system.append(Atom(sym_inner, position=(cx, cy, cz)))
    else:
        # Single atom gas proxy
        system.append(Atom(sym_outer, position=(cx, cy, cz)))

    # Setup LAMMPS full atom style requirements (Mol-ID and Charges)
    nmols = np.ones(len(system), dtype=int)
    nmols[n_graphene_atoms:] = 2 # Mark new atoms as Molecule 2
    
    system.set_array('mol-id', nmols)
    system.set_initial_charges(np.zeros(len(system)))

    # -------------------------------------------------------------------
    # SMART FIX: Lock the Type Order globally using specorder
    # Type 1 = C, Type 2 = H, Type 3 = F, Type 4 = Cl
    # -------------------------------------------------------------------
    write(args.output, system, format='lammps-data', atom_style='full', specorder=['C', 'H', 'F', 'Cl'])
    
    # Ultimate safeguard: Force '4 atom types' in the file header
    with open(args.output, 'r') as f:
        lines = f.readlines()
    with open(args.output, 'w') as f:
        for line in lines:
            if "atom types" in line:
                f.write("4 atom types\n")
            else:
                f.write(line)
    
    print(f"Success: Added {args.gas} using robust proxy mapping. Saved to {args.output}")

if __name__ == "__main__":
    main()
