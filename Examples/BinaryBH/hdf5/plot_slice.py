import numpy as np
#import matplotlib.pyplot as plt
import h5py

f = h5py.File('BinaryBHChk_000000.3d.hdf5', 'r')
base_items = list(f.items())
print("Items in the base directory:", base_items)
