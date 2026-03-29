import numpy as np
from ase.io import read, write
from ase.build import make_supercell

def create_4nm_graphene(input_cif):
    atoms = read(input_cif)
    single_layer = atoms[[atom.index for atom in atoms if atom.position[2] < 2.0]]

    P = [[2, 1, 0], [0, 1, 0], [0, 0, 1]]
    orthogonal_unit = make_supercell(single_layer, P)
    
    lengths = orthogonal_unit.cell.lengths()
    orthogonal_unit.set_cell([lengths[0], lengths[1], 100.0], scale_atoms=True)
    
    target_length = 40.0
    nx = int(np.round(target_length / lengths[0]))
    ny = int(np.round(target_length / lengths[1]))
    
    supercell = orthogonal_unit.repeat((nx, ny, 1))
    
    # ---> FIX: Force the graphene layer exactly to the middle (Z = 50.0) <---
    z_shift = 50.0 - supercell.positions[:, 2].mean()
    supercell.positions[:, 2] += z_shift
    
    write("graphene_4nm_base.xyz", supercell)

if __name__ == "__main__":
    create_4nm_graphene("Graphite.cif")
