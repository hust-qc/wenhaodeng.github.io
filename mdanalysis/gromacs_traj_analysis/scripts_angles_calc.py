import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis import rms

# set font
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']

u = mda.Universe('step3_input.gro', 'step7_npt2_production.xtc')

# please select three atoms using index name
atom1 = u.select_atoms(f"index 31717")  # O5'
atom2 = u.select_atoms(f"index 31718")  # PA, this is the center of the angle you selected
atom3 = u.select_atoms(f"index 32200")  # O1

# calculated angles
angles = []
for ts in u.trajectory:
    # this is a vector 
    vec1 = (atom2.positions - atom1.positions).flatten()  # one line of the angle
    vec2 = (atom2.positions - atom3.positions).flatten()  # the other line of the angle
    
    # math 
    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)
    if norm1 == 0 or norm2 == 0:
        angles.append(np.nan)  
        continue
    cos_angle = np.dot(vec1, vec2) / (norm1 * norm2)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    
    # degree
    angle = np.arccos(cos_angle) * (180.0 / np.pi)
    angles.append(angle)
angles = np.array(angles)

#draw figure 
plt.figure(figsize=(10, 6))
plt.hist(valid_angles, bins=25, color='purple', edgecolor='black', alpha=0.7, label='Angle Distribution')
plt.xlabel('Angle (degrees)', fontsize=16)
plt.ylabel('Frequency', fontsize=16)
# plt.title('Histogram of Angle between Three Atoms', fontsize=16)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)
plt.tick_params(axis='both', labelsize=12)

plt.xlim(120, 180)  # please modify the range
plt.ylim(0, 200)   # please modify the range

plt.savefig('angle_histogram_chaina_o5pa_o1.png', dpi=600)
plt.close()
