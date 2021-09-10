import numpy as np
#import matplotlib.pyplot as plt
import h5py

f = h5py.File('FixedBoostLL_plt000290.3d.hdf5', 'r')
base_items = list(f.items())
print("Items in the base directory:", base_items)
G1 = f.get("Chombo_global")
G1_items = list(G1.items())
print("Items in G1:", G1_items)
