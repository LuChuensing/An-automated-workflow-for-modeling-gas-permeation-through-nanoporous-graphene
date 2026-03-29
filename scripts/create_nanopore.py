import numpy as np
from ase.io import read, write
from ase.neighborlist import neighbor_list
from ase import Atom

def extract_hollow_and_passivate_hole_only(input_file, n_limit, output_file):
    atoms = read(input_file)
    atoms.pbc = True 
    
   
    ii, jj = neighbor_list('ij', atoms, cutoff=1.8, self_interaction=False)
    adj = [set() for _ in range(len(atoms))]
    neighbors_old = [[] for _ in range(len(atoms))]
    for i, j in zip(ii, jj):
        adj[i].add(j)
        neighbors_old[i].append(j)

    
    all_rings = []
    for i in range(len(atoms)):
        for n1 in adj[i]:
            for n2 in adj[n1] - {i}:
                for n3 in adj[n2] - {n1}:
                    for n4 in adj[n3] - {n2}:
                        for n5 in adj[n4] - {n3}:
                            if i in adj[n5]:
                                all_rings.append(tuple(sorted([i, n1, n2, n3, n4, n5])))
    unique_rings = [tuple(r) for r in set(all_rings)]
    num_rings = len(unique_rings)
    
    pos = atoms.get_positions()
    ring_centers = np.array([np.mean(pos[list(r)], axis=0) for r in unique_rings])

    
    sheet_center = np.mean(pos, axis=0)
    start_idx = int(np.argmin(np.linalg.norm(ring_centers - sheet_center, axis=1)))
    
    # (Dual Graph)
    ring_adj = [[] for _ in range(num_rings)]
    for i in range(num_rings):
        set_i = set(unique_rings[i])
        for j in range(i + 1, num_rings):
            if len(set_i.intersection(unique_rings[j])) == 2:
                ring_adj[i].append(j)
                ring_adj[j].append(i)

    
    ordered_rings = []
    visited = set([start_idx])
    queue = [start_idx]
    
    while queue and len(ordered_rings) < n_limit:
        curr = queue.pop(0)
        ordered_rings.append(curr)
        
        neighbors = [n for n in ring_adj[curr] if n not in visited]
        if neighbors:
            rel_vecs = ring_centers[neighbors] - ring_centers[start_idx]
            angles = np.arctan2(rel_vecs[:, 1], rel_vecs[:, 0])
            sorted_neighbors = [n for _, n in sorted(zip(angles, neighbors))]
            for n in sorted_neighbors:
                visited.add(n)
                queue.append(n)

    ordered_rings = ordered_rings[:n_limit]

    
    atom_to_rings = [set() for _ in range(len(atoms))]
    for r_idx, r_atoms in enumerate(unique_rings):
        for a_idx in r_atoms:
            atom_to_rings[a_idx].add(r_idx)

    to_delete = []
    ordered_set = set(ordered_rings)
    for a_idx, r_ids in enumerate(atom_to_rings):
        if len(r_ids) >= 2 and r_ids.issubset(ordered_set):
            to_delete.append(a_idx)

    print(f"RESULT")
    print(f"Target Area S: {n_limit}")
    
    if len(to_delete) > 0:
        
        to_delete_set = set(to_delete)
        edge_tags = np.zeros(len(atoms), dtype=int)
        
        for i in range(len(atoms)):
            if i not in to_delete_set:
               
                if any(n in to_delete_set for n in neighbors_old[i]):
                    edge_tags[i] = 1 # 打上标记 1
                    
        
        atoms.set_tags(edge_tags)

        
        del atoms[to_delete]
        print(f"Step 1: Successfully removed {len(to_delete)} Carbon atoms.")
        
        
        
        ii_new, jj_new = neighbor_list('ij', atoms, cutoff=1.8, self_interaction=False)
        neighbors_new = [[] for _ in range(len(atoms))]
        for i, j in zip(ii_new, jj_new):
            neighbors_new[i].append(j)
            
        tags_remaining = atoms.get_tags()
        h_atoms_to_add = []
        C_H_LENGTH = 1.09 
        
        for i, neighs in enumerate(neighbors_new):
            
            if len(neighs) == 2 and tags_remaining[i] == 1:
                v1 = atoms.get_distance(i, neighs[0], mic=True, vector=True)
                v2 = atoms.get_distance(i, neighs[1], mic=True, vector=True)
                
                v1_dir = v1 / np.linalg.norm(v1)
                v2_dir = v2 / np.linalg.norm(v2)
                
                
                v_out = -(v1_dir + v2_dir)
                norm_out = np.linalg.norm(v_out)
                
                if norm_out > 0.1:
                    v_out_dir = v_out / norm_out
                    h_pos = atoms.positions[i] + v_out_dir * C_H_LENGTH
                    
                    h_atoms_to_add.append(Atom('H', position=h_pos))
                    
        
        if h_atoms_to_add:
            for h in h_atoms_to_add:
                atoms.append(h)
            print(f"Step 2: Successfully added {len(h_atoms_to_add)} H atoms to passivate the pore. External boundaries are untouched.")
        
        write(output_file, atoms)
        print(f"File saved to: {output_file}")
    else:
        print("Warning: No internal atoms were identified for removal.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input base graphene xyz")
    parser.add_argument("-n", "--size", type=int, required=True, help="Nanopore area S")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    args = parser.parse_args()

    extract_hollow_and_passivate_hole_only(args.input, args.size, args.output)
