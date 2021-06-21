# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# output data for setup
M = 1.0
symmetry = 8
# make the plot
fig = plt.figure()

# volume integral dataset
data1 = np.loadtxt("VolumeIntegrals.dat")
timedata = data1[:,0]
dM = symmetry*data1[:,1] - symmetry*data1[0,1]
Source = symmetry*data1[:,3]

# flux dataset integration
data1 = np.loadtxt("SurfaceIntegrals.dat")
timedata = data1[:,0]
NetEoFlux = data1[:,1]
f_flux = interpolate.interp1d(timedata, NetEoFlux, kind = 'cubic')
f_source = interpolate.interp1d(timedata, Source, kind = 'cubic')
finetimedata = np.linspace(timedata[0], timedata[np.size(timedata)-1], 10000)
dt = finetimedata[1] - finetimedata[0]
FEodt = np.zeros_like(finetimedata)
Sourcedt = np.zeros_like(finetimedata)
for i, t in enumerate(finetimedata) :
    FEodt[i] += FEodt[i-1] + f_flux(t) * dt
    Sourcedt[i] += Sourcedt[i-1] + f_source(t) * dt

plt.figure(figsize=(8,6))
plt.plot(timedata, dM, '-', lw = 2.0, label=r'$\int (Q-Q_0) dV$')
plt.plot(finetimedata, Sourcedt, '--', lw = 1.5, label=r'$\int \int \mathcal{S} dV dt$')
plt.plot(finetimedata, FEodt, ':', lw = 2.0, label=r'$\int \int F dS dt$')
# make the plot look nice
plt.rc('axes', titlesize=16)
plt.xlabel("simulation time")
plt.ylabel("energy (LL)")
plt.xlim(0, 82)
plt.ylim(-0.022, 0.022)
plt.rc('legend', fontsize=14)
plt.legend(loc=0)
plt.grid()

# save as png image
filename = "EvsT_LL.png"
plt.tight_layout()
plt.savefig(filename)
