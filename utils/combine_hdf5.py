import h5py
import numpy as np
import sys
import os

def copy_dataset(src_ds, dest_group, name):
    if name not in dest_group:
        dest_group.create_dataset(name, data=src_ds[()], dtype=src_ds.dtype)

def copy_group(src_group, dest_group):
    for name, item in src_group.items():
        if isinstance(item, h5py.Dataset):
            copy_dataset(item, dest_group, name)
        elif isinstance(item, h5py.Group):
            if name not in dest_group:
                new_group = dest_group.create_group(name)
                copy_group(item, new_group)
            else:
                copy_group(item, dest_group[name])

def merge_hdf5_files(file1, file2, output_file):
    with h5py.File(file1, 'r') as f1, h5py.File(file2, 'r') as f2, h5py.File(output_file, 'w') as fout:
        # Copy everything from file1 (takes priority)
        for name in f1:
            f1.copy(name, fout)

        # Copy from file2 only if not already present
        for name in f2:
            if name not in fout:
                f2.copy(name, fout)
            else:
                copy_group(f2[name], fout[name])

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merge_h5.py <output_file> <input_file_1> <input_file_2>")
        sys.exit(1)

    output_file = sys.argv[1]
    input_file1 = sys.argv[2]
    input_file2 = sys.argv[3]

    if not os.path.isfile(input_file1) or not os.path.isfile(input_file2):
        print("Error: One or both input files do not exist.")
        sys.exit(1)

    merge_hdf5_files(input_file1, input_file2, output_file)
    print(f"Merged HDF5 file created at: {output_file}")


hf = h5py.File(str(sys.argv[2]), 'r')
