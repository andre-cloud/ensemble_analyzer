import matplotlib.pyplot as plt
import numpy as np 


uv_ref = np.loadtxt('/Users/andrea/Desktop/scratch/vit_c/ECD_ref_norm.xy')
# p3 = np.loadtxt('/Users/andrea/Desktop/scratch/vit_c/UV_p3_comp.xy')
p4 = np.loadtxt('/Users/andrea/Desktop/scratch/vit_c/ECD_p4_comp.xy')
p4[:,1] /= np.max(np.abs(p4[:,1]))


plt.plot(uv_ref[:, 0], uv_ref[:, 1])
# plt.plot(p3[:, 0], p3[:, 1])
plt.plot(p4[:, 0], p4[:, 1])
plt.show()