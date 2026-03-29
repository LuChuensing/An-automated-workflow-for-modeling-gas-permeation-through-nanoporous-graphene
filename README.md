# GraphenePore-ZScan
An automated Python-based workflow for modeling graphene nanopores and calculating gas permeation energy barriers via Z-direction potential energy scanning in LAMMPS.
## Installation
Ensure you have Python and LAMMPS installed.The workflow requires the ASE.
```
pip install ase
conda install -c conda-forge lammps
```
## Workflow Overview
This pipeline automates the entire simulation setup in four main stages:  
1.System Generation: Converts a graphite CIF file into a passivated monolayer graphene sheet with a central nanopore.  
2.Molecular Assembly: Automatically identifies the pore center and inserts a gas molecule (Ar, CH4, N2, H2, or CO2) at a specified distance.  
3.Forcefield Mapping: Applies Lorentz-Berthelot mixing rules to generate precise non-bonded interaction parameters between the gas and the graphene/hydrogen atoms.  
4.Z-Scan Simulation: Generates and executes a LAMMPS script that moves the gas molecule through the pore, recording the interaction energy at each step to create a potential energy profile.  
## File Structure & Functions
`master.py`: The central controller. Orchestrates the execution of all sub-scripts in the correct order.  
`build_monolayer.py`: Extracts a single layer from a bulk CIF/data file and expands it into a supercell.  
`create_nanopore.py`: Deletes carbon atoms within a radius and passivates the broken bonds with Hydrogen.  
`add_mol.py`: Inserts the gas molecule and assigns specific Atom Types for LAMMPS.
`make_input.py`: produce input file for LAMMPS.
## Usage
Place your Graphite.cif in the project directory and run the master script:
```
python master.py
```
The results will be saved in a directory
