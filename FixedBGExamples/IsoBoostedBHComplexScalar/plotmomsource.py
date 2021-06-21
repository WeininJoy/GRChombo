# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# output data for setup
M = 1.0
symmetry = 4
# make the plot
fig = plt.figure()

# volume integral dataset
data1 = np.loadtxt("VolumeIntegrals.dat")
timedata = data1[:,0]
dM = symmetry*data1[:,3] - symmetry*data1[0,3]
Source = symmetry*data1[:,4]

# flux dataset integration
data1 = np.loadtxt("SurfaceIntegrals.dat")
timedata = data1[:,0]
NetEoFlux = data1[:,6]
NetEiFlux = data1[:,3]
fo_flux = interpolate.interp1d(timedata, NetEoFlux, kind = 'cubic')
fi_flux = interpolate.interp1d(timedata, NetEiFlux, kind = 'cubic')
f_source = interpolate.interp1d(timedata, Source, kind = 'cubic')
finetimedata = np.linspace(timedata[0], timedata[np.size(timedata)-1], 10000)
dt = finetimedata[1] - finetimedata[0]
FEodt = np.zeros_like(finetimedata)
FEidt = np.zeros_like(finetimedata)
Sourcedt = np.zeros_like(finetimedata)
for i, t in enumerate(finetimedata) :
    FEodt[i] += FEodt[i-1] + fo_flux(t) * dt
    FEidt[i] += FEidt[i-1] + fi_flux(t) * dt
    Sourcedt[i] += Sourcedt[i-1] + f_source(t) * dt

plt.figure(figsize=(8,6))
plt.plot(finetimedata, FEidt - FEodt + Sourcedt, '--', lw = 2.5, label=r'$\int (\int(F_i - F_o) dS + \int \mathcal{S} dV) dt$')
plt.plot(timedata, dM, '--', lw = 1.5, label=r'$\int (Q-Q_0) dV$')
plt.plot(finetimedata, Sourcedt, ':', lw = 2.0, label=r'$\int \int \mathcal{S} dV dt$')
plt.plot(finetimedata, FEodt, ':', lw = 2.0, label=r'$\int \int F_o dS dt$')
plt.plot(finetimedata, FEidt, ':', lw = 2.0, label=r'$\int \int F_i dS dt$')

# make the plot look nice
plt.rc('axes', titlesize=16)
plt.xlabel("simulation time")
plt.ylabel("momentum")
#plt.xlim(0, 82)
#plt.ylim(-0.022, 0.022)
plt.rc('legend', fontsize=14)
plt.legend(loc=0)
plt.grid()

# save as png image
filename = "MvsT_DF.png"
plt.tight_layout()
plt.savefig(filename)
