import ROOT
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as st

def get_every_fourth_entry_starting_at_first(arr):
    """
    Returns every fourth element in a list, starting from the first element.

    Args:
        arr: The input list.

    Returns:
        A new list containing every fourth element of the input list, 
        starting from the first, or an empty list if the input list is empty.
    """

    if not arr:
        return []

    return arr[0::4]  # Slice from the 1st element (index 0) with a step of 4



# Example usage:
# my_array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
# result = get_every_fourth_entry_starting_at_first(my_array)
# print(result)  # Output: [1, 5, 9, 13]

# empty_array = []
# result_empty = get_every_fourth_entry_starting_at_first(empty_array)
# print(result_empty) # Output: []

# short_array = [1, 2, 3]
# result_short = get_every_fourth_entry_starting_at_first(short_array)
# print(result_short) # Output: [1]

# short_array_2 = [1, 2, 3, 4]
# result_short_2 = get_every_fourth_entry_starting_at_first(short_array_2)
# print(result_short_2) # Output: [1]



def plot_top_mass(root_file="2LSS_prediction_v01.root", tree_name="reco", branch_name="top_m_SPANET", whichtop="any",
                minPlot = 0, maxPlot= 500, ):
    """
    Plots the specified branch from a ROOT file within a given range.

    Args:
        root_file (str): Path to the ROOT file.
        tree_name (str): Name of the TTree.
        branch_name (str): Name of the branch to plot.
    """

    try:
        # Open the ROOT file
        root_file_obj = ROOT.TFile.Open(root_file, "READ")
        if not root_file_obj or root_file_obj.IsZombie():
            raise Exception(f"Error opening file: {root_file}")

        # Get the TTree
        tree = root_file_obj.Get(tree_name)
        if not tree:
            raise Exception(f"Error getting tree: {tree_name}")

        for i in range(min(5, tree.GetEntries())): #Check up to 5 entries
            tree.GetEntry(i)
            print(getattr(tree, "top_m_SPANET")) #Pri
        # Get the branch
        branch = tree.GetBranch(branch_name)
        if not branch:
            raise Exception(f"Error getting branch: {branch_name}")
        print("YO")
        print(branch)

        # Efficiently read the branch into a numpy array (recommended)
        n_entries = tree.GetEntries()
        top_mass_values = np.zeros(4*n_entries, dtype=np.float32)  # Use float32 for memory efficiency
        truth_top_mass_values = np.zeros(4*n_entries, dtype=np.float32)  # Use float32 for memory efficiency
        parton_top_isHadronic_values = np.zeros(4*n_entries, dtype=np.float32)  # Use float32 for memory efficiency
        print(top_mass_values)
        print(n_entries)
        for i in range(n_entries):
            truth_top_isHad_mass_values_vector = getattr(tree, "parton_top_isHadronic")
            print(len(truth_top_isHad_mass_values_vector))
            print("TEST")
            print(i)
            tree.GetEntry(i)
            print(tree)
            # assign_prob_values_vector = getattr(tree, "parton_top_m") 
            top_mass_values_vector = getattr(tree, branch_name) # More flexible way to access branch value
            j = 0
            for  topMass in top_mass_values_vector:
                top_mass_values[i*4+j] = topMass #/1000
                j = j +1
            truth_top_mass_values_vector = getattr(tree, "parton_top_m") 
            truth_top_isHad_mass_values_vector = getattr(tree, "parton_top_isHadronic") 
            j = 0
            for  topMass in truth_top_mass_values_vector:
                parton_top_isHadronic_values[i*4+j] = truth_top_isHad_mass_values_vector[j] 
                if truth_top_isHad_mass_values_vector[j] == 1:
                    truth_top_mass_values[i*4+j] = topMass/1000
                else:
                    truth_top_mass_values[i*4+j]  = -0.001
                j = j +1


        len_of_isHad = []
        for i in range(n_entries):
            tree.GetEntry(i)
            truth_top_isHad_mass_values_vector = getattr(tree, "parton_top_isHadronic")
            print(len(truth_top_isHad_mass_values_vector))
            print(truth_top_isHad_mass_values_vector)
            len_of_isHad.append(sum(truth_top_isHad_mass_values_vector))
        # plt.hist(len_of_isHad, bins=6, range=(-1, 5),  histtype='step', color='blue') #histtype='step' for a line plot
        # plt.savefig(f"comparison_len_of_isHad_2LSS_wLep_{whichtop}_v03.pdf") #Save the plot        
                
                
                
                                    
        print(len(top_mass_values))
        print(len(parton_top_isHadronic_values))
        # print(top_mass_values)
        # OR, if you don't have enough memory for very large files, iterate and append:
        #top_mass_values = []
        #for entry in tree:
        #    top_mass_values.append(getattr(entry, branch_name))
        #top_mass_values = np.array(top_mass_values)

        print()
        # Define the plotting range (0 to 1.4^6)
        min_mass = minPlot
        max_mass =maxPlot  # or 1.4e6 if you meant 1.4 * 10^6
        # Apply the range cut to the numpy array
        # top_mass_values_cut = top_mass_values[1: :4]# + top_mass_values[1::4]
        
        t1 = top_mass_values[0: :4]
        t2 = top_mass_values[1: :4]
        t3 = top_mass_values[2: :4]
        t4 = top_mass_values[3: :4]
        
        print(len(len_of_isHad))
        print(len_of_isHad)
        print(len(t4))
        len_of_isHad = np.array(len_of_isHad)
        # t1 = t1[len_of_isHad==2]
        # t2 = t2[len_of_isHad==2]
        # t3 = t3[len_of_isHad==2]
        # t4 = t4[len_of_isHad==2]

        t1 = t1[len_of_isHad!=2]
        t2 = t2[len_of_isHad!=2]
        t3 = t3[len_of_isHad!=2]
        t4 = t4[len_of_isHad!=2]         
        print(len(t4))   
        # exit()
##
        
        if (whichtop == "any"):
            top_mass_values_cut = np.append(t1,t2)
            top_mass_values_cut = np.append(top_mass_values_cut,t3)
            top_mass_values_cut = np.append(top_mass_values_cut,t4)

        elif (whichtop== "hadronic"):
            top_mass_values_cut = np.append(t1,t2)
            
        elif (whichtop =="leptonic"):
            top_mass_values_cut = np.append(t3,t4)
            
        
        # top_mass_values_cut = top_mass_values_cut[(top_mass_values_cut > min_mass)]        
        # truth_top_mass_values = truth_top_mass_values[(truth_top_mass_values > min_mass) & (truth_top_mass_values <= max_mass)]
        
        
        print("Len of array being plotted: ", len(top_mass_values_cut))        
        
        # Create the histogram using matplotlib (more flexible)
        top_mass_value = 172.76
        W_mass_value = 80.37
        mbl_limit = np.sqrt((top_mass_value*top_mass_value -W_mass_value*W_mass_value)*(W_mass_value*W_mass_value))/W_mass_value
        mbl_limit = round(mbl_limit,2)
        

        plt.hist(top_mass_values_cut, bins=50, range=(min_mass, max_mass),  histtype='step', color='blue', label=branch_name) #histtype='step' for a line plot
        if (branch_name == 'top_m_SPANET'):
            top_mass_value = 172.76
            plt.axvline(x=top_mass_value *1e3, color='red', linestyle='--', label=f"Top Mass ({top_mass_value} GeV)")  # Added line
            # plt.hist(truth_top_mass_values, bins=50, range=(min_mass, max_mass),  histtype='step', color='green', label="Truth") #histtype='step' for a line plot
        if (branch_name == 'W_m_SPANET'):
            W_mass_value = 80.37
            plt.axvline(x=W_mass_value *1e3, color='red', linestyle='--', label=f"W Mass ({W_mass_value} GeV)")  # Added line
        if (branch_name == 'mbl_SPANET'):
            print(mbl_limit)
            # W_mass_value = mbl_limit
            plt.axvline(x=mbl_limit*1e3, color='red', linestyle='--', label=f"MBL LIMIT = {mbl_limit} GeV)")  # Added line
        # plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))



        print("plotted")
        plt.title(f"Distribution of {branch_name}")
        plt.xlabel(branch_name)
        plt.ylabel("Frequency")
        plt.legend()
        plt.grid(True)
        # plt.yscale('log') #Optional: Use a logarithmic y-axis if needed
        plt.tight_layout() # Adjusts subplot params so that subplots fit in to the figure area.
        plt.savefig(f"not2LSS_comparison_{branch_name}_2LSS_wLep_{whichtop}_v03.pdf") #Save the plot
        plt.show()


    except Exception as e:
        print(f"An error occurred: {e}")

    finally:
        if root_file_obj:
            root_file_obj.Close()  # Important: Close the ROOT file

# Example usage:
# plot_top_mass() #Uses default values, or you can specify the file, tree, and branch names:
# plot_top_mass( branch_name="mbl_SPANET")
# plot_top_mass(root_file = "testtest.root",branch_name="W_m_SPANET")
# plot_top_mass(root_file = "testtest.root")
# plot_top_mass(root_file = "testtest.root",branch_name="W_m_SPANET")
# plot_top_mass(branch_name="W_m_SPANET") #2LSS_prediction_v01
# plot_top_mass(root_file = "2LSS_Hadronic_prediction_v07.root")
root_file = "/hpcfs/groups/phoenix-hpc-coepp/atlas/mjgreen/four-top/SPANet/2LSS_wLep_v03.root"
# top mass
# plot_top_mass(root_file =root_file,minPlot = 0, maxPlot = 500e3,whichtop = "hadronic")


# W mass
# plot_top_mass(root_file = root_file,branch_name="W_m_SPANET")
plot_top_mass(root_file =root_file,branch_name="W_m_SPANET",minPlot = -2, maxPlot = 500e3,whichtop = "hadronic")

# mbl
# plot_top_mass(root_file =root_file,branch_name="mbl_SPANET", minPlot = 0, maxPlot = 400e3,whichtop = "leptonic")

# Assignment Prob
# plot_top_mass(root_file =root_file,branch_name="top_assign_prob_SPANET",whichtop = "any", minPlot = 0, maxPlot = 1)



# b_recondex
# plot_top_mass(root_file =root_file,branch_name="b_recoj_index", minPlot = -2, maxPlot = 10,whichtop = "leptonic")

# b_recondex
# plot_top_mass(root_file =root_file,branch_name="wd1_recoj_index", minPlot = -2, maxPlot = 10,whichtop = "parton_top_m")


# plot_top_mass(root_file =root_file,branch_name="parton_top_m", minPlot = -2, maxPlot = 500e3,whichtop = "leptonic")