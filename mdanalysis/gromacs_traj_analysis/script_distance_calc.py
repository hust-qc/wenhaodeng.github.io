import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis import rms

# 设置全局字体为Arial
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']

u = mda.Universe('step3_input.gro', 'step7_npt2_production.xtc')

#check the index of the selected atoms
atom1 = u.select_atoms(f"index 32200")  #o1 hydroxide
atom2 = u.select_atoms(f"index 31718")   #pa of dATP
print(atom1)
print(atom2)

if len(atom1) != 1 or len(atom2) != 1:
    raise ValueError("One or both atoms not found. Please check the index numbers.")

# add each distance into a list
distances = []  
for ts in u.trajectory:
    distance = np.linalg.norm(atom1.positions - atom2.positions)
    distances.append(distance)

# to be an array
distances = np.array(distances)

#draw figure using plt
plt.figure(figsize=(10, 6))
plt.hist(distances, bins=25, color='red', edgecolor='black', range=(2.6, 3.4), alpha=0.7, label='Distance Distribution')  # please change range if you want
plt.xlabel('Distance (Å)', fontsize=16)
plt.ylabel('Frequency', fontsize=16)
# plt.title('Histogram of Atom Distance')
plt.grid(True, alpha=0.3)
plt.legend(fontsize=12)

plt.tick_params(axis='both', labelsize=12)

# set limitations
plt.xlim(2.6, 3.4)  # please change this according to your select bond
plt.ylim(0, 200)   # please change this according to the frequency

# save image
plt.savefig('atom_distance_histogram_chaina_o1pa_mn.png', dpi=600)
plt.close()
