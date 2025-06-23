import h5py
import sys
import ROOT
import numpy as np
import glob
import matplotlib.pyplot as plt
from train_test_splitter import  filter_and_save_h5,count_events,filter_and_save_h5_random
import math
import os


def delta_phi(phi1, phi2):
    """Calculate the difference in phi between two angles, wrapped to [-π, π]."""
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2 * math.pi
    while dphi < -math.pi:
        dphi += 2 * math.pi
    return dphi

def delta_r(eta1, phi1, eta2, phi2):
    """Calculate delta R between two particles with (eta, phi) coordinates."""
    deta = eta1 - eta2
    dphi = delta_phi(phi1, phi2)
    return math.sqrt(deta**2 + dphi**2)


def resolve_conflict(index_list_a, deltaR_list_a, index_list_b, deltaR_list_b, verbose=False):
    for i, idx in enumerate(index_list_a):
        if idx != -1 and idx in index_list_b:
            j = index_list_b.index(idx)
            if deltaR_list_a[i] <= deltaR_list_b[j]:
                index_list_b[j] = -1
            else:
                if verbose:
                    print(deltaR_list_a[i])
                index_list_a[i] = -1

def resolve_internal_duplicates(index_list, deltaR_list):
    for i, idx in enumerate(index_list):
        if idx == -1:
            continue
        for j in range(i + 1, len(index_list)):
            if index_list[j] == idx:
                if deltaR_list[i] <= deltaR_list[j]:
                    index_list[j] = -1
                else:
                    index_list[i] = -1
                    break  # No need to keep checking once i is invalidated

def match_partons_to_jets(
    parton_eta, parton_phi, pdg_ids, jet_eta, jet_phi, 
    reco_index_new, deltaR_new, deltaR_tolerance,  
    wd1_el_index_new  = None, wd1_mu_index_new=None, nElectrons=None, reco_index_original=None, max_deltaR_holder=None, 
    skip_neutrinos=False, label="", DEBUG = False
):
    for i in range(len(parton_eta)):
        if skip_neutrinos and abs(pdg_ids[i]) in (12, 14,16):
            continue

        min_deltaR = .4
        min_index = -1
        eta1, phi1 = parton_eta[i], parton_phi[i]

        for j in range(len(jet_eta)):
            if jet_eta[j] == 0:
                continue
            dR = delta_r(eta1, phi1, jet_eta[j], jet_phi[j])
            if dR < min_deltaR:
                min_deltaR = dR
                min_index = j

        if min_deltaR < deltaR_tolerance:
            reco_index_new[i] = min_index
            deltaR_new[i] = min_deltaR

        if DEBUG:
            print(f"{label} minimum_deltaR: {min_deltaR}, index: {min_index}", end="")
            if pdg_ids is not None:
                print(f", pdgId: {pdg_ids[i]}")
            else:
                print()

        # Check for lepton mapping to 16 or 17
        if pdg_ids is not None and abs(pdg_ids[i]) in (11, 13,15):  # electron or muon or tau 
            if min_index not in (16, 17):
                # print(f"⚠️  Warning: {label} lepton (pdgId {pdg_ids[i]}) matched to suspicious jet index {min_index}")
                reco_index_new[i] = -1
                continue
            else:
                if nElectrons == 0:
                    if min_index ==16:
                        wd1_el_index_new[i] = 0
                    else:
                        wd1_el_index_new[i] = 1
                if nElectrons == 1:
                    if min_index ==16:
                        wd1_el_index_new[i] = 0
                    else:
                        wd1_mu_index_new[i] = 0
                if nElectrons == 2:
                    if min_index ==16:
                        wd1_mu_index_new[i] = 0
                    else:
                        wd1_mu_index_new[i] = 1                        

                    
                

        if pdg_ids is not None and abs(pdg_ids[i]) not in (11, 13,15):  # electron or muon or tau
            if min_index in (16, 17):
                # print(f"⚠️  Warning: {label} jet (pdgId {pdg_ids[i]}) matched to suspicious lepton index {min_index}")
                reco_index_new[i] = -1
                continue           
        if (label=="b" and min_index in (16, 17)):
                # print(f"⚠️  Warning: {label} jet  matched to suspicious lepton index {min_index}")
                reco_index_new[i] = -1
                # continue    
            
                
        # Update max deltaR if needed
        if (
            max_deltaR_holder is not None and 
            reco_index_original is not None and
            min_index == reco_index_original[i] and
            min_deltaR > max_deltaR_holder[0]
        ):
            max_deltaR_holder[0] = min_deltaR


def pt_eta_phi_e_to_px_py_pz_e(pt, eta, phi, E):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    return np.array([E, px, py, pz])


def calculate_invariant_mass(*quarks):
    """
    Calculates invariant mass given 1 to 3 quarks in (pt, eta, phi, E) format.

    Args:
        quarks: variable number of quark 4-vectors, each a tuple (pt, eta, phi, E)

    Returns:
        Invariant mass (float)
    """
    if len(quarks) == 0 or len(quarks) > 3:
        raise ValueError("Must provide between 1 and 3 quark four-vectors.")

    total_vec = np.zeros(4)  # [E, px, py, pz]

    for q in quarks:
        total_vec += pt_eta_phi_e_to_px_py_pz_e(*q)

    mass_squared = total_vec[0]**2 - np.dot(total_vec[1:], total_vec[1:])
    return np.sqrt(mass_squared) if mass_squared > 0 else 0.0

import matplotlib.pyplot as plt
from collections import Counter

def plot_pie_chart_from_array(data,plotDirectory,prefix):
    """
    Plots a pie chart showing the percentage of 1s, 2s, and 3s in the input array.
    
    Parameters:
        data (list of int): List containing only 1s, 2s, and 3s.
    """
    # Count occurrences
    counts = Counter(data)

    # Only include 1, 2, and 3 in case others are present
    filtered_counts = {k: counts[k] for k in (0,1, 2, 3) if k in counts}

    if not filtered_counts:
        print("No 1s, 2s, or 3s in the input data.")
        return

    # Prepare data
    labels = [str(k) for k in filtered_counts.keys()]
    sizes = list(filtered_counts.values())

    # Plot
    plt.figure(figsize=(6, 6))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
    plt.title('Percentage of 0s,1s, 2s, and 3s')
    plt.axis('equal')
    plt.savefig(plotDirectory +"/"+prefix +"_nJets.pdf") #Save the plot

    
def pt_eta_phi_to_pxyz(pt, eta, phi):
    """
    Convert pt, eta, phi arrays to px, py, pz arrays.

    Parameters:
        pt (np.ndarray): Transverse momentum array
        eta (np.ndarray): Pseudorapidity array
        phi (np.ndarray): Azimuthal angle array

    Returns:
        tuple: (px, py, pz) as separate np.ndarrays
    """
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    return px, py, pz