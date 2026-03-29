import os
import sys
import subprocess
import shutil

def run_cmd(cmd, cwd=None):
    """Helper function to execute shell commands securely."""
    print(f"🚀 Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        print(f"❌ Error during execution of: {' '.join(cmd)}")
        sys.exit(1)

def main():
    print("=== Graphene Nanopore Workflow Manager ===")
    
    if not os.path.exists("Graphite.cif"):
        print("❌ Error: 'Graphite.cif' not found in the current directory.")
        sys.exit(1)
    
    print("\n--- Step 1: Building Monolayer ---")
    run_cmd([sys.executable, "build_monolayer.py"])
    base_xyz = "graphene_4nm_base.xyz"

    try:
        s_val = int(input("\nPlease enter the desired nanopore area S (e.g., 10): ").strip())
    except ValueError:
        print("❌ Invalid input. Must be an integer.")
        sys.exit(1)
        
    pore_dir = f"Nanopore{s_val}"
    os.makedirs(pore_dir, exist_ok=True)
    pore_xyz = os.path.join(pore_dir, f"graphene_N{s_val}_perfect_hole.xyz")
    
    print("\n--- Step 2: Creating Nanopore ---")
    run_cmd([sys.executable, "create_nanopore.py", "-i", base_xyz, "-n", str(s_val), "-o", pore_xyz])

    valid_gases = ['H2', 'N2', 'CH4', 'Ar', 'CO2']
    gas_name = input(f"\nPlease enter the gas molecule name ({'/'.join(valid_gases)}): ").strip()
    if gas_name not in valid_gases:
        print(f"❌ Invalid gas name. Must be one of {valid_gases}.")
        sys.exit(1)
        
    gas_dir = os.path.join(pore_dir, gas_name)
    os.makedirs(gas_dir, exist_ok=True)
    data_file = os.path.join(gas_dir, f"system_{gas_name}.data")
    
    print("\n--- Step 3: Adding Gas Molecule ---")
    run_cmd([sys.executable, "add_mol.py", "-i", pore_xyz, "-g", gas_name, "-o", data_file])

    print("\n--- Step 4: Generating LAMMPS Pull Script ---")
    run_cmd([sys.executable, "make_input.py", data_file])

    print("\n--- Step 5: Running LAMMPS Simulation ---")
    # Automatically copy the AIREBO potential file into the execution directory
    if os.path.exists("CH.airebo"):
        shutil.copy("CH.airebo", os.path.join(gas_dir, "CH.airebo"))
    else:
        print("⚠️ Warning: 'CH.airebo' not found in root dir. LAMMPS might crash if it needs it.")
        
    print("⏳ This may take a few moments...")
    run_cmd(["lmp", "-in", "pull.in"], cwd=gas_dir)

    print("\n✅ Workflow & Simulation Complete!")
    print(f"➡️  Check '{gas_dir}/barrier_{gas_name}.txt' for your results!")

if __name__ == "__main__":
    main()
