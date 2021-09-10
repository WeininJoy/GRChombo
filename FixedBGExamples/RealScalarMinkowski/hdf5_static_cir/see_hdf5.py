import numpy as np
import h5py

f = h5py.File('FixedBoostLL_plt000000.3d.hdf5', 'r')
level_0 = f['level_0']
data = level_0['data:datatype=0']
print(data[:])
