import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis import rms

# set Arial font for all text in the generated figure
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']

# import the setup model and the trajectory of 100 ns production, 
u = mda.Universe('step3_input.gro', 'step7_npt2_production.xtc') 

# select backbone as reference
protein = u.select_atoms('protein and backbone')

# crystal structure as the base
ref = mda.Universe('step3_input.gro')
ref_protein = ref.select_atoms('protein and backbone')

# calculate RMSD value
rmsd_analysis = rms.RMSD(protein, ref_protein, select='backbone')
rmsd_analysis.run()

# 100 ns here we used, please assure how long the production runs 
time = rmsd_analysis.times/1000  #please print (rmsd_analysis.times) to check the time it runs and modify this value when necessary!
rmsd = rmsd_analysis.rmsd[:, 2]   

# calculate average RMSD
average_rmsd = np.mean(rmsd)
print(f"Average RMSD: {average_rmsd:.3f} Å")

# plot
plt.figure(figsize=(10, 6))
plt.plot(time, rmsd, 'b-', label='RMSD')
plt.xlabel('Time (ns)', fontsize=16)
plt.ylabel('RMSD (Å)', fontsize=16)
#plt.title('RMSD over Time')
plt.grid(True)
plt.legend(fontsize=14)
plt.tick_params(axis='both', labelsize=14)

# save image
plt.savefig('rmsd_plot_mn.png', dpi=600)
plt.close()
