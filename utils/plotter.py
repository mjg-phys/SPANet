import h5py
import sys
import numpy as np
from array import array
import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as st
 

def plotterFunction(data= [],
                       labels = [],
                       nameOfVariables = "topMass_hadronic",
                       xlim=None, 
                       ylim=None,
                       plotter_directory = "plots/",
                       plotWMass = False,
                       plotTMass = False,
                       plot_mbl_limit = False,
                       doLog = True,
                       colourArray = None
                       ):

    # Create the histogram using matplotlib (more flexible)
    top_mass_value = 172.76
    W_mass_value = 80.37
    
    xmin = xlim[0]
    xmax= xlim[1]
    mbl_limit = np.sqrt((top_mass_value*top_mass_value -W_mass_value*W_mass_value)*(W_mass_value*W_mass_value))/W_mass_value
    mbl_limit = round(mbl_limit,2)
    print(data)
    
    for thisFile,thisLabel,thisColour in zip(data,labels,colourArray):
        print(thisLabel)
        print(len(thisFile))
        print(min(thisFile))
        plt.hist(thisFile, bins=50, range=(xmin, xmax),  histtype='step', color=thisColour, label=thisLabel) #histtype='step' for a line plot
    
    if (plotTMass):
        top_mass_value = 172.76
        plt.axvline(x=top_mass_value *1e3, color='red', linestyle='--', label=f"Top Mass ({top_mass_value} GeV)")  # Added line
        # plt.hist(truth_top_mass_values, bins=50, range=(min_mass, max_mass),  histtype='step', color='green', label="Truth") #histtype='step' for a line plot
    if (plotWMass):
        W_mass_value = 80.37
        plt.axvline(x=W_mass_value *1e3, color='red', linestyle='--', label=f"W Mass ({W_mass_value} GeV)")  # Added line
    if (plot_mbl_limit):
        plt.axvline(x=mbl_limit*1e3, color='red', linestyle='--', label=f"MBL LIMIT = {mbl_limit} GeV)")  # Added line

    plt.title(f"Distribution of {nameOfVariables}")
    plt.xlabel(nameOfVariables)
    plt.ylabel("Frequency")
    plt.legend()
    plt.grid(True)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    # plt.yscale('log') #Optional: Use a logarithmic y-axis if needed
    plt.tight_layout() # Adjusts subplot params so that subplots fit in to the figure area.
    plt.savefig( plotter_directory+nameOfVariables+"_v4" +".pdf") #Save the plot
    plt.show()
    
    
    


def plot_2d_histogram(data, title='2D Histogram', xlabel='Predicted', ylabel='True',
                      bins=50, cmap='viridis', xlim=None, ylim=None,plotter_directory="plots/",
                      plotTMass=True,plotWMass=True):
    """
    Plots a 2D histogram from two arrays of equal length.

    Parameters:
        data (array-like): A list or array of two arrays (x, y), both of same length.
        title (str): Title of the plot.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        bins (int or tuple): Number of bins or a tuple (xbins, ybins).
        cmap (str): Colormap used for the plot.
        xlim (tuple): Tuple (xmin, xmax) to set x-axis limits.
        ylim (tuple): Tuple (ymin, ymax) to set y-axis limits.
    """
    if len(data) != 2:
        raise ValueError("Input 'data' must be a list or array of two arrays.")

    x, y = np.array(data[0]), np.array(data[1])

    if x.shape != y.shape:
        raise ValueError("Both arrays must be of the same shape.")

    # Create 2D histogram
    plt.figure(figsize=(8, 6))
    hist = plt.hist2d(x, y, bins=bins, cmap=cmap)

    # Apply axis limits
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)

    # Add colorbar and labels
    plt.colorbar(label='Counts')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if (plotTMass):
        top_mass_value = 172.76
        plt.axvline(x=top_mass_value *1e3, color='red', linestyle='--', label=f"Top Mass ({top_mass_value} GeV)")  # Added line
        plt.axhline(y=top_mass_value *1e3, color='red', linestyle='--', label=f"Top Mass ({top_mass_value} GeV)")  # Added line
        
    if (plotWMass):
        top_mass_value = 80.37
        plt.axvline(x=top_mass_value *1e3, color='red', linestyle='--', label=f"Top Mass ({top_mass_value} GeV)")  # Added line
        plt.axhline(y=top_mass_value *1e3, color='red', linestyle='--', label=f"Top Mass ({top_mass_value} GeV)")  # Added line


    # plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig( plotter_directory+title +"2D.pdf") #Save the plot
    
    
import sys
import h5py

# Skip the script name, take all other arguments as filenames
input_filenames = sys.argv[1:]

# Open each file in r+ mode and store in a list
hf_array = [h5py.File(fname, 'r+') for fname in input_filenames]


# hf = h5py.File(str(sys.argv[1]), 'r+')


# truth_topMass

plottingVector = []

####################################################################################

# plottingTypes = [
#     # ["REGRESSIONS/t1/tMass","REGRESSIONS/t2/tMass" ],
#             # ["REGRESSIONS/t1/tMass_SPANET","REGRESSIONS/t2/tMass_SPANET" ]
#                 ["TRUTH/EVENT/t3_nu_cos_phi"],
#     ["REGRESSIONS/EVENT/t3_nu_cos_phi" ],
# ]
# nameOfVariable = "nu_cos_phi"
# labels = ["Truth", "SPANet"]
# colour_array = ["red", "green"]
# # x_max = 3
# # x_min = -3
# xlim = [-1.5, 1.5]

# nameOfVariable = "t2_tMass"



# nameOfVariable = "mu_adjustedTruth"
# plottingTypes = [
#     # ["REGRESSIONS/t1/tMass","REGRESSIONS/t2/tMass" ],
#             # ["REGRESSIONS/t1/PARTICLE/tMass","REGRESSIONS/t2/PARTICLE/tMass"],
#             ["REGRESSIONS/t1/wMass","REGRESSIONS/t2/wMass"],
#     # ["REGRESSIONS/EVENT/"+nameOfVariable],
#     ["REGRESSIONS/t1/wMass_SPANET","REGRESSIONS/t2/wMass_SPANET"],
#     # ["REGRESSIONS/t1/tMass_SPANET","REGRESSIONS/t2/tMass_SPANET"],
#     # ["REGRESSIONS/t1/tMass_SPANET","REGRESSIONS/t2/tMass_SPANET"],
#     ["TRUTH/EVENT/t1_wMass","TRUTH/EVENT/t2_wMass"],

# ]


# nameOfVariable = "t3_nu_pz"
# plottingTypes = [
#     ["TRUTH/EVENT/t3_nu_pz" ],
#             ["REGRESSIONS/EVENT/t3_nu_pz"],
# ]
# xlim = [-3, 3]


# nameOfVariable = "tHadronic_mass"
# plottingTypes = [
#      ["TRUTH/EVENT/t1_tMass","TRUTH/EVENT/t2_tMass"],
#     ["REGRESSIONS/t1/tMass_SPANET","REGRESSIONS/t2/tMass_SPANET"],
#     ["REGRESSIONS/t1/PARTICLE/tMass","REGRESSIONS/t2/PARTICLE/tMass"],
#     ["REGRESSIONS/EVENT/t1_tMass","REGRESSIONS/EVENT/t2_tMass"],

# ]
# xlim = [0, 400000]


# nameOfVariable = "wHadronic_mass"
# plottingTypes = [
#      ["TRUTH/EVENT/t1_wMass","TRUTH/EVENT/t2_wMass"],
#     ["REGRESSIONS/t1/wMass_SPANET","REGRESSIONS/t2/wMass_SPANET"],
#     ["REGRESSIONS/t1/wMass","REGRESSIONS/t2/wMass"],
#     ["REGRESSIONS/EVENT/t1_wMass","REGRESSIONS/EVENT/t2_wMass"],

# ]
# xlim = [0, 400000]




nameOfVariable = "mbl"
plottingTypes = [
    ["REGRESSIONS/t3/mbl","REGRESSIONS/t4/mbl"],
    ["REGRESSIONS/t3/mbl_SPANET","REGRESSIONS/t4/mbl_SPANET"],
    ["TRUTH/EVENT/t3_mbl","TRUTH/EVENT/t3_mbl"],
    ["REGRESSIONS/EVENT/t3_mbl","REGRESSIONS/EVENT/t3_mbl"],

]
xlim = [-1000, 400000]



labels = plottingTypes

fileIndices = [0,0,0,0,0,0,0,0,0,0]
# labels = ["Original Truth", 
#           "No Regression", 
#         #   "+ regression",
#         #   "+ non-zero scale regression",
#           "Post Cleaning"]
colour_array = ["red", "green", "purple","blue","brown"]
# x_max = 3
# x_min = -3
# xlim = [0, 400000]



# plottingTypes = [
#     ["REGRESSIONS/t1/wMass_SPANET","REGRESSIONS/t2/wMass_SPANET" ],
#         ["REGRESSIONS/t1/wMass","REGRESSIONS/t2/wMass" ]

# ]
# nameOfVariable = "wMass_hadronic"

# labels = ["Truth", "SPANet"]
# colour_array = ["red", "green"]
# x_max = 400000


# plottingTypes = [
#     ["REGRESSIONS/t3/mbl","REGRESSIONS/t4/mbl" ],
#     ["REGRESSIONS/t3/mbl_SPANET","REGRESSIONS/t4/mbl_SPANET" ],
# ]
# nameOfVariable = "mbl_hadronic"

# labels = ["Truth", "SPANet"]
# colour_array = ["red", "green"]
# x_max = 400000



fullTypes = []
for thisLoc,thisIndex in zip(plottingTypes,fileIndices):
    arrayToPlot = []
    for thisType in thisLoc: 
        print(thisType)
        arrayToPlot.append(np.array(hf_array[thisIndex].get(thisType)))
        # print(arrayToPlot)
    arrayToPlot = np.concatenate(arrayToPlot)
    fullTypes.append(arrayToPlot)   
    
     



plotterFunction(fullTypes,labels,xlim=xlim,nameOfVariables=nameOfVariable, colourArray= colour_array)

# plot_2d_histogram(fullTypes,
#                   title=nameOfVariable,
#                   xlabel='Predicted Value',
#                   ylabel='True Value',
#                   bins=300,
#                   xlim =[-1,1],
#                   ylim =[-1,1])