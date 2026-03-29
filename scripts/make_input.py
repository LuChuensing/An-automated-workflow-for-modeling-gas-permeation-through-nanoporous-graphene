import argparse
import sys
import os
import math

def generate_pull_script(data_file):
    if not os.path.exists(data_file):
        print(f"Error: {data_file} not found.")
        sys.exit(1)

    gas_name = None
    for g in ['H2', 'N2', 'CO2', 'CH4', 'Ar']:
        if g in data_file:
            gas_name = g
            break
            
    if not gas_name:
        print("Error: Could not determine gas type from filename.")
        sys.exit(1)


    K_TO_EV = 8.617333e-5 

    # Parameters: [ (mass_outer, eps_outer_K, sig_outer), (mass_inner, eps_inner_K, sig_inner) ]
    gas_params = {
        'H2':  [ (1.008,  0.0,   0.0),    (0.0001, 36.7,  2.958) ],
        'N2':  [ (14.007, 36.0,  3.31),   (0.0001, 0.0,   0.0)   ],
        'CO2': [ (15.999, 79.0,  3.05),   (12.011, 27.0,  2.80)  ],
        'CH4': [ (16.043, 148.0, 3.73) ],
        'Ar':  [ (39.948, 125.7, 3.345) ]
    }

    params = gas_params[gas_name]
    is_3_site = len(params) == 2

    # Graphene-C and Edge-H LJ Parameters (converted from K to eV)
    mass_C, eps_C, sig_C = 12.011, 28.0  * K_TO_EV, 3.40
    mass_H, eps_H, sig_H = 1.008,  25.45 * K_TO_EV, 2.36

    # Extract Gas Parameters in eV
    mass_3 = params[0][0]
    eps_3 = params[0][1] * K_TO_EV
    sig_3 = params[0][2]

    if is_3_site:
        mass_4 = params[1][0]
        eps_4 = params[1][1] * K_TO_EV
        sig_4 = params[1][2]
        gas_grp = "3 4"
    else:
        # Dummy parameters for single-site gases (Type 4 won't interact)
        mass_4 = 1.0
        eps_4 = 0.0
        sig_4 = 0.0
        gas_grp = "3"

    # =====================================================================
    # Manual Lorentz-Berthelot Mixing Rules for Hybrid Cross Interactions
    # eps_ij = sqrt(eps_i * eps_j)  (Now in eV)
    # sig_ij = (sig_i + sig_j) / 2
    # =====================================================================
    eps_13 = math.sqrt(eps_C * eps_3)
    sig_13 = (sig_C + sig_3) / 2.0
    eps_23 = math.sqrt(eps_H * eps_3)
    sig_23 = (sig_H + sig_3) / 2.0

    eps_14 = math.sqrt(eps_C * eps_4)
    sig_14 = (sig_C + sig_4) / 2.0
    eps_24 = math.sqrt(eps_H * eps_4)
    sig_24 = (sig_H + sig_4) / 2.0

    # Self-cross term for the gas (Type 3 - Type 4)
    eps_34 = math.sqrt(eps_3 * eps_4)
    sig_34 = (sig_3 + sig_4) / 2.0

    # =====================================================================
    # Constructing the LAMMPS Script
    # =====================================================================
    script = f"""# Advanced Z-Scan Script with AIREBO Relaxation for {gas_name}

units           metal
atom_style      full
boundary        p p f

# 1. HYBRID FORCE FIELD DEFINITION
pair_style      hybrid airebo 3.0 lj/cut 12.0

read_data       {os.path.basename(data_file)}

# 2. MASSES
mass            1 {mass_C}
mass            2 {mass_H}
mass            3 {mass_3}
mass            4 {mass_4}

# 3. FORCE FIELD ASSIGNMENTS
# AIREBO applied strictly to C (Type 1) and H (Type 2)
pair_coeff      * * airebo CH.airebo C H NULL NULL

# Gas-Gas pure LJ interactions (Energy in eV)
pair_coeff      3 3 lj/cut {eps_3:.8f} {sig_3:.6f}
pair_coeff      4 4 lj/cut {eps_4:.8f} {sig_4:.6f}
pair_coeff      3 4 lj/cut {eps_34:.8f} {sig_34:.6f}

# Manual Cross-Interactions (Graphene <-> Gas)
pair_coeff      1 3 lj/cut {eps_13:.8f} {sig_13:.6f}
pair_coeff      2 3 lj/cut {eps_23:.8f} {sig_23:.6f}
pair_coeff      1 4 lj/cut {eps_14:.8f} {sig_14:.6f}
pair_coeff      2 4 lj/cut {eps_24:.8f} {sig_24:.6f}

# 4. GROUPS
group           gas type {gas_grp}
group           graphene type 1 2

# Exclude gas-gas internal LJ interactions if you treat it as rigid later
neigh_modify    exclude group gas gas

# =====================================================================
# PHASE 1: RELAXATION OF GRAPHENE PORE
# =====================================================================
fix             freeze_gas gas setforce 0.0 0.0 0.0

thermo          10
thermo_style    custom step pe ke etotal press
print           "--- Starting Graphene Pore Relaxation ---"
minimize        1.0e-8 1.0e-6 10000 100000

unfix           freeze_gas
print           "--- Relaxation Complete ---"

# =====================================================================
# PHASE 2: Z-SCAN PULLING
# =====================================================================
fix             freeze_gra graphene setforce 0.0 0.0 0.0

compute         inter_pe gas group/group graphene
print           "Z_Position Interaction_PE_eV" file barrier_{gas_name}.txt

variable        z_start equal 65.0
variable        z_end equal 35.0
variable        dz equal -0.1
variable        nsteps equal abs((v_z_end-v_z_start)/v_dz)

print           "--- Starting Z-Scan Pulling ---"

variable        i loop ${{nsteps}}
label           scan_loop
    displace_atoms  gas move 0.0 0.0 ${{dz}}
    fix             lock_gas gas setforce 0.0 0.0 0.0
    run             0 
    variable        current_z equal xcm(gas,z)
    variable        gas_wall_pe equal c_inter_pe
    print           "${{current_z}} ${{gas_wall_pe}}" append barrier_{gas_name}.txt
    unfix           lock_gas
next            i
jump            SELF scan_loop

print           "--- Simulation Complete ---"
"""

    out_dir = os.path.dirname(data_file)
    pull_path = os.path.join(out_dir, "pull.in") if out_dir else "pull.in"
    with open(pull_path, "w", encoding="utf-8") as f:
        f.write(script)
    print(f"Success: Generated advanced {pull_path} with METAL units and AIREBO relaxation.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file")
    args = parser.parse_args()
    generate_pull_script(args.data_file)
