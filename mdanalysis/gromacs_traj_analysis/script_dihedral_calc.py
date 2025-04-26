import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis import rms

# set font of text
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']

u = mda.Universe('step3_input.gro', 'step7_npt2_production.xtc')

atom1 = u.select_atoms(f"index 31714")  # C5'
atom2 = u.select_atoms(f"index 31717")  # O5'
atom3 = u.select_atoms(f"index 31718")  # PA
atom4 = u.select_atoms(f"index 32200")  # O1

# calculate dihehydral
dihedrals = []
for ts in u.trajectory:
    # from postion to a array [3,]
    p1 = atom1.positions.flatten()  #atom1.positions is a NumPy two dimention array
    p2 = atom2.positions.flatten()
    p3 = atom3.positions.flatten()
    p4 = atom4.positions.flatten()
    
    # the normal vectors of two planes
    # the first one
    vec1 = p2 - p1
    vec2 = p3 - p2
    normal1 = np.cross(vec1, vec2)
    # the second one
    vec3 = p3 - p2
    vec4 = p4 - p3
    normal2 = np.cross(vec3, vec4)
    # check
    norm1 = np.linalg.norm(normal1)
    norm2 = np.linalg.norm(normal2)
    if norm1 == 0 or norm2 == 0:
        dihedrals.append(np.nan)  
        continue
    
    # math
    cos_dihedral = np.dot(normal1, normal2) / (norm1 * norm2)
    cos_dihedral = np.clip(cos_dihedral, -1.0, 1.0)
    dihedral = np.arccos(cos_dihedral) * (180.0 / np.pi)
    triple_product = np.dot(normal1, np.cross(normal2, p3 - p2))
    if triple_product < 0:
        dihedral = -dihedral
    
    dihedrals.append(dihedral)

dihedrals = np.array(dihedrals)

valid_dihedrals = dihedrals[~np.isnan(dihedrals)]
if len(valid_dihedrals) > 0:
    average_dihedral = np.mean(valid_dihedrals)
    print(f"Average Dihedral Angle between planes: {average_dihedral:.3f} degrees")
else:
    print("No valid dihedral angles calculated.")

# draw figure
plt.figure(figsize=(10, 6))
plt.hist(valid_dihedrals, bins=50, color='blue', edgecolor='black', alpha=0.7, label='Dihedral Angle Distribution')
plt.xlabel('Dihedral Angle (degrees)', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.title('Histogram of Dihedral Angle', fontsize=16)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)
plt.tick_params(axis='both', labelsize=12)

#set range
plt.xlim(-180, 180)
plt.ylim(0, 800)

# save image
plt.savefig('dihedral_histogram_c5p_o5p_pa_o1.png', dpi=600)
plt.close()
