import h5py
import numpy as np
from sklearn.model_selection import train_test_split

# Load the .h5 file
# file_path = '/hpcfs/groups/phoenix-hpc-coepp/atlas/mjgreen/four-top/run/new_samples/test.h5'  # Replace with your file path
# with h5py.File(file_path, 'r') as file:
#     # Check the keys in the .h5 file (datasets inside)
#     print("Keys in the HDF5 file:", list(file.keys()))

#     # Access the dataset safely
#     data = file['INPUTS']  # Don't use [:] yet, just get the dataset reference

#     # if isinstance(data, h5py.Dataset):
#     #     print("Shape of INPUTS:", data.shape)  # Get dataset shape
#     #     print("Data type of INPUTS:", data.dtype)  # Get dataset dtype
#     # else:
#     #     print("INPUTS is not a dataset but a group.")
        
#     if isinstance(data, h5py.Group):
#         print("Keys inside INPUTS:", list(data.keys()))
        
#     leptons = data['Lepton']
#     if isinstance(leptons, h5py.Group):
#         print("Keys inside INPUTS:", list(leptons.keys()))
    
#     pt_lep = leptons["pt"]
#     print(pt_lep)
    
    
import h5py
import numpy as np

def filter_and_save_h5(input_file, train_file, val_file):
    with h5py.File(input_file, 'r') as f:
        assigned_jets = f['INPUTS/Event/assigned_jets'][:]

        # Create boolean mask
        train_mask = (assigned_jets >= 4).squeeze()
        val_mask = ~train_mask  # Complement mask

        # Open new files to save filtered data
        with h5py.File(train_file, 'w') as train_f, h5py.File(val_file, 'w') as val_f:
            
            def copy_filtered_data(source_group, dest_group, mask):
                """ Recursively copy datasets while applying the mask. """
                for key, item in source_group.items():
                    print(f"Processing key: {key}")  # Debug print
                    if isinstance(item, h5py.Group):
                        print(f"Creating group: {key}")  # Debug print
                        new_group = dest_group.create_group(key)
                        copy_filtered_data(item, new_group, mask)  # Recursive call
                    elif isinstance(item, h5py.Dataset):
                        print(f"Dataset: {key}, Shape: {item.shape}, Mask Shape: {mask.shape}")  # Debug print
                        if item.shape[0] == mask.shape[0]:  # Ensure shape compatibility
                            print(f"Filtering dataset: {key}")  # Debug print
                            filtered_data = np.array(item)[mask]  # Convert to numpy array and apply mask
                            dest_group.create_dataset(key, data=filtered_data)
                        else:
                            print(f"Copying dataset without filtering: {key}")  # Debug print
                            dest_group.create_dataset(key, data=item[:])  # Copy without filtering if shape mismatch
            
            # Copy the structure while filtering data
            copy_filtered_data(f, train_f, train_mask)
            copy_filtered_data(f, val_f, val_mask)

def count_events(file_path):
    with h5py.File(file_path, 'r') as f:
        num_events = f['INPUTS/Event/assigned_jets'].shape[0]  # First dimension = number of events
        print(f['INPUTS/Event/assigned_jets'][:])
        return num_events

def filter_and_save_h5_random(input_file, train_file, val_file, train_fraction=0.8):
    with h5py.File(input_file, 'r') as f:
        assigned_jets = f['INPUTS/Event/assigned_jets'][:]

        # Generate shuffled indices
        num_events = assigned_jets.shape[0]
        indices = np.arange(num_events)
        # np.random.shuffle(indices)

        # Split indices based on fraction
        split_idx = int(num_events * train_fraction)
        train_indices = indices[:split_idx]
        val_indices = indices[split_idx:]

        # Open new files to save filtered data
        with h5py.File(train_file, 'w') as train_f, h5py.File(val_file, 'w') as val_f:
            
            def copy_filtered_data(source_group, dest_group, indices):
                """ Recursively copy datasets while applying shuffled indices. """
                for key, item in source_group.items():
                    print(f"Processing key: {key}")  # Debug print
                    if isinstance(item, h5py.Group):
                        print(f"Creating group: {key}")  # Debug print
                        new_group = dest_group.create_group(key)
                        copy_filtered_data(item, new_group, indices)  # Recursive call
                    elif isinstance(item, h5py.Dataset):
                        print(f"Dataset: {key}, Shape: {item.shape}")  # Debug print
                        if item.shape[0] == num_events:  # Ensure shape compatibility
                            print(f"Filtering dataset: {key}")  # Debug print
                            filtered_data = np.array(item)[indices]  # Apply shuffled indices
                            dest_group.create_dataset(key, data=filtered_data)
                        else:
                            print(f"Copying dataset without filtering: {key}")  # Debug print
                            dest_group.create_dataset(key, data=item[:])  # Copy without filtering if shape mismatch
            
            # Copy the structure while filtering data
            copy_filtered_data(f, train_f, train_indices)
            copy_filtered_data(f, val_f, val_indices)



def checkFile(file_path):
    with h5py.File(file_path, 'r') as f:
        num_events = f['INPUTS/Event/assigned_jets'].shape[0]  # First dimension = number of events
        bjets = f['INPUTS/Event/assigned_bjets'][:]  # First dimension = number of events
        jets = f['INPUTS/Event/assigned_jets'][:]  # First dimension = number of events
        qjets = f['INPUTS/Event/assigned_q_jets'][:]  # First dimension = number of events
        bjects = f['INPUTS/Event/assigned_objects'][:]  # First dimension = number of events
        njets = f['INPUTS/Event/nJets'][:]  # First dimension = number of events
        jet_eta = f['INPUTS/Jet/eta'][:]  # First dimension = number of events
        jet_pt = f['INPUTS/Jet/pt'][:]  # First dimension = number of events

        print(file_path)
        print("num_events: ", num_events)
        print("njets: ", njets)
        print("jet_pt: ", jet_pt)
        numJets = 0
        numJets_v = []
        for i in range(len(jet_eta)):
            numJets = 0
            for j in range(len(jet_eta[i])):
                if jet_eta[i][j] != 0:
                    numJets = numJets+1
            numJets_v.append(numJets)
        print(numJets_v)
        # print("bjets: ", bjets)
        # print("qjets", qjets)
        # print("jets:", jets)
        # print("objects", bjects)
        # print(f['INPUTS/Event/assigned_jets'][:])
        print(len(njets))
        print(len(numJets_v))
        return num_events
    
    
    
# filter_and_save_h5('old_all_samples_v2.h5', 'train.h5', 'val.h5')
# test_events = count_events("train.h5")
# val_events = count_events("val.h5")

# # print(f"Number of events in test.h5: {test_events}")
# # print(f"Number of events in val.h5: {val_events}")
# print(f"Ratio train.h5: {test_events/(test_events+val_events)}")
# # print(f"Number of events in val.h5: {val_events}")

# # Example usage
# filter_and_save_h5_random('old_all_samples_v2.h5', 'train_random.h5', 'val_random.h5', train_fraction=test_events/(test_events+val_events))

# test_events = count_events("train_random.h5")
# val_events = count_events("val_random.h5")

# print(f"Number of events in test.h5: {test_events}")
# print(f"Number of events in val.h5: {val_events}")
# print(f"Ratio train.h5: {test_events/(test_events+val_events)}")
# print(f"Number of events in val.h5: {val_events}")


checkFile("train_geq4jets.h5")

checkFile("val_geq4jets.h5")

checkFile("train_random.h5")

checkFile("train_original.h5")


checkFile("testtest.h5")