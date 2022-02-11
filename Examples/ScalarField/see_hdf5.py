import h5py

filename = "hdf5/ScalarFieldp_000000.3d.hdf5"

with h5py.File(filename, "r") as f:
    # List all groups
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[1]

    # Get the data
    print( list(f[a_group_key].keys()) )
    sub_group_key = list(f[a_group_key].keys())[0]

    data = f[a_group_key][sub_group_key]
    print(data)
