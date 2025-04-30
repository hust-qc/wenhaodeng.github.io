#this script could be used for calculations of rmsd/ angle/ distance/dihedral angle
#if you want to change the color for the generated figures, please read this script and change plt code part.
import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis import rms
import os
from MDAnalysis.analysis.dihedrals import Dihedral

# font Arial
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']

# Function to check if a file exists
def check_file_exists(filename):
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File '{filename}' not found. Please check the file name or path.")

# Interactive input for topology and trajectory files
print("Enter the topology file (e.g., mode_input.gro，mode_input.pdb) and press Enter:")
topology_file = input().strip()
check_file_exists(topology_file)

print("Enter the trajectory file (e.g., npt_production.xtc,npt_production.dcd) and press Enter:")
trajectory_file = input().strip()
check_file_exists(trajectory_file)

# Load the universe with user-provided files
try:
    u = mda.Universe(topology_file, trajectory_file)
except Exception as e:
    raise RuntimeError(f"Failed to load universe with topology '{topology_file}' and trajectory '{trajectory_file}': {str(e)}")

# Select atoms for RMSD calculation (protein backbone)
protein = u.select_atoms('protein and backbone')

# Set reference structure (using the topology file)
try:
    ref = mda.Universe(topology_file)
    ref_protein = ref.select_atoms('protein and backbone')
except Exception as e:
    raise RuntimeError(f"Failed to load reference structure from '{topology_file}': {str(e)}")

# Function to validate and get atom index
def get_atom_index(prompt, universe):
    while True:
        try:
            index = int(input(prompt).strip())
            # Check if the index is valid by attempting to select the atom
            atom = universe.select_atoms(f"index {index}")
            if not atom:
                print(f"No atom found with index {index}. Please enter a valid index.")
                continue
            return atom
        except ValueError:
            print("Please enter a valid integer for the atom index.")
# Prompt for analysis type
while True:
    analysis_shot = input("Enter analysis type (distance, angle, dihedral) and press Enter: ").strip().lower()
    if analysis_shot in ['distance', 'angle', 'dihedral']:
        break
    print("Invalid analysis type. Please enter 'distance', 'angle', or 'dihedral'.")

# Determine number of atoms needed
num_atoms = {'distance': 2, 'angle': 3, 'dihedral': 4}[analysis_shot]

# Get atom indices
atoms = []
for i in range(num_atoms):
    prompt = f"Enter index for atom {i+1} and press Enter: "
    atom = get_atom_index(prompt, u)
    atoms.append(atom)

#RMSD analysis
rmsd_analysis = rms.RMSD(protein, ref_protein, select='backbone')
rmsd_analysis.run()

# Extract time and RMSD data
time = rmsd_analysis.times / 1000  # Convert to ns, PLEASE change this value based on the timestep you set
rmsd = rmsd_analysis.rmsd[:, 2]    # RMSD values

# Calculate average RMSD
average_rmsd = np.mean(rmsd)
print(f"Average RMSD: {average_rmsd:.3f} Å")

# Plot RMSD curve
plt.figure(figsize=(10, 6))
plt.plot(time, rmsd, 'b-', label='RMSD')
plt.xlabel('Time (ns)', fontsize=18)
plt.ylabel('RMSD (Å)', fontsize=18)
plt.grid(True)
plt.legend(fontsize=16)
plt.tick_params(axis='both', labelsize=16)

# Save the figure RMSD
output_file = 'rmsd_plot_traj.png'
plt.savefig(output_file)
plt.close()
print(f"RMSD plot saved as '{output_file}'")

# Analyze based on analysis type
results = []
if analysis_shot == 'distance':
    # Calculate distance between two atoms
    for ts in u.trajectory:
        dist = mda.lib.distances.calc_bonds(atoms[0].positions, atoms[1].positions, box=u.dimensions)
        results.append(dist[0])
    results = np.array(results)  # In angstroms
    ylabel = 'Distance (Å)'
    title = 'Distance Histogram'
    output_file = 'distance_histogram.png'

elif analysis_shot == 'angle':
    # Calculate angle between three atoms
    for ts in u.trajectory:
        angle = mda.lib.distances.calc_angles(
            atoms[0].positions, atoms[1].positions, atoms[2].positions, box=u.dimensions
        )
        results.append(np.degrees(angle[0]))  # Convert to degrees
    results = np.array(results)
    ylabel = 'Angle'
    title = 'Angle Histogram'
    output_file = 'angle_histogram.png'

else:  # dihedral
    # Calculate dihedral angle between four atoms
    dihedral_atoms =atoms[0] + atoms[1] + atoms[2] + atoms[3]  #universe
    dihedral = Dihedral([dihedral_atoms]).run()
    results = dihedral.results.angles.flatten()  # In degrees
    ylabel = 'Dihedral Angle'
    title = 'Dihedral Histogram'
    output_file = 'dihedral_histogram.png'

# Create figure histogram
plt.figure(figsize=(10, 6))
plt.hist(results, bins=25, color='mediumseagreen', alpha=0.7, edgecolor='black')  #please change the color if you want
plt.xlabel(ylabel, fontsize=18)
plt.ylabel('Frequency', fontsize=18)
plt.title(title, fontsize=16)
plt.grid(True, alpha=0.3)
plt.tick_params(axis='both', labelsize=16)

# Save the figure
plt.savefig(output_file, dpi=600)
plt.close()
print(f"Histogram saved as '{output_file}'")
