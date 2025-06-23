import h5py
import sys
import ROOT
import numpy as np
import glob
import matplotlib.pyplot as plt
from train_test_splitter import  filter_and_save_h5,count_events,filter_and_save_h5_random
from helper_functions import delta_phi, delta_r, resolve_conflict
from helper_functions import resolve_internal_duplicates, match_partons_to_jets
from helper_functions import calculate_invariant_mass,plot_pie_chart_from_array
from helper_functions import pt_eta_phi_to_pxyz
maxjets = 16
import math
import os

chain = ROOT.TChain("reco")

maximum_deltaR_bparton = 0
maximum_deltaR_Wd1 = 0
maximum_deltaR_Wd2 = 0
nEvent = 100000000
doDeltaRMatching = True 
do_fuzzy_matching = True
doDEBUG = False
addTruthInformation = True
prefix = "newMatching_v8"

# Only keep the branches that you need for the preprocessing so it doesnt take forever. 
with open("/hpcfs/groups/phoenix-hpc-coepp/atlas/mjgreen/four-top/SPANet/treesToKeep.txt") as f:
    branches_to_keep = [line.strip() for line in f if line.strip()]


# Get the ROOT files
folder_path = sys.argv[2]  # Make sure this is a directory, not a file
root_files = glob.glob(f"{folder_path}/*.root")
for root_file in root_files:
    chain.Add(root_file)

print(f"Added {len(root_files)} files to the chain.")


# Uncomment this to get all the branches
# for branch in chain.GetListOfBranches():
#     print(branch.GetName())
# exit()

# Disable all branches and enable only desired branches
chain.SetBranchStatus("*", 0)
for branch in branches_to_keep:
    chain.SetBranchStatus(branch, 1)

cut_2LSS = "(pass_SSee_NOSYS||pass_SSem_NOSYS||pass_SSmm_NOSYS)"
cut_3L = "(pass_eee_ZVeto_NOSYS||pass_eem_ZVeto_NOSYS||pass_emm_ZVeto_NOSYS||pass_mmm_ZVeto_NOSYS)"
cut = None


do2LSS = True
do3L = False

if do2LSS:
    cut = cut_2LSS
else if do3L:
    cut = cut_3L
else:
    print("Please select either 2LSS or 3L")
    exit()

# Define the cut and do it
tree = chain.CopyTree(cut,"", nEvent)

name = str(sys.argv[1])

if '.h5' not in name:
    name += '.h5'


running_count= []

emptylist = [[0]*maxjets]*tree.GetEntries()

# Fill data lists
mask, eta, etag, btag, mtag, phi, pt, e, q = [], [], [],[], [], [], [], [], []
jet_mask, jet_eta, jet_etag, jet_btag, jet_mtag, jet_phi, jet_pt, jet_e, jet_q = [], [], [],[], [], [], [], [], []
jet_mask, jet_eta, jet_etag, jet_btag, jet_mtag, jet_phi, jet_pt, jet_e, jet_q = [], [], [],[], [], [], [], [], []

lepton_mask, lepton_eta, lepton_etag, lepton_btag, lepton_mtag, lepton_phi, lepton_pt, lepton_e, lepton_q = [], [], [],[], [], [], [], [], []
nJets_ii = []

b2, l2, b3, l3, b4, l4 = [], [], [], [], [], []
b1, q11, q12, b2, q21, q22 = [], [], [], [], [], []
kk = 0

number_matched_vector = []
met_met = []
met_sin_phi = []
met_cos_phi = []

leptonic_W_d1_pdg_id= []
leptonic_W_d2_pdg_id= []
W_leptonic_lep3_truth_m,W_leptonic_lep3_truth_pt,W_leptonic_lep3_truth_eta,W_leptonic_lep3_truth_phi= [],[],[],[]
W_leptonic_nu3_truth_m,W_leptonic_nu3_truth_pt,W_leptonic_nu3_truth_eta,W_leptonic_nu3_truth_phi= [],[],[],[]

W_leptonic_lep4_truth_m,W_leptonic_lep4_truth_pt,W_leptonic_lep4_truth_eta,W_leptonic_lep4_truth_phi= [],[],[],[]
W_leptonic_nu4_truth_m,W_leptonic_nu4_truth_pt,W_leptonic_nu4_truth_eta,W_leptonic_nu4_truth_phi= [],[],[],[]

W_hadronic_q11_truth_m,W_hadronic_q11_truth_pt,W_hadronic_q11_truth_eta,W_hadronic_q11_truth_phi= [],[],[],[]
W_hadronic_q12_truth_m,W_hadronic_q12_truth_pt,W_hadronic_q12_truth_eta, W_hadronic_q12_truth_phi= [],[],[],[]

W_hadronic_q21_truth_m,W_hadronic_q21_truth_pt,W_hadronic_q21_truth_eta,W_hadronic_q21_truth_phi= [],[],[],[]
W_hadronic_q22_truth_m,W_hadronic_q22_truth_pt,W_hadronic_q22_truth_eta, W_hadronic_q22_truth_phi= [],[],[],[]

t3_nu_eta = []
t3_nu_cos_phi = []
t3_nu_sin_phi = []

t4_nu_eta = []
t4_nu_cos_phi = []
t4_nu_sin_phi = []

leptonic_W_d2_pdg_id= []
print("Starting the prepocessing")



for event in tree:        
    if len(event.parton_top_isHadronic) != 4:
        continue
    if sum(event.parton_top_isHadronic) != 2:
        continue
    kk = kk+ 1
    if kk % 1000 == 0:
        print(f"Processed {kk} events")    
    lepton_maski, lepton_etai, lepton_etagi, lepton_mtagi, lepton_btagi, lepton_phii, lepton_pti, lepton_ei, lepton_qi = [], [], [],[], [], [], [], [], []
    
    jet_maski, jet_etai, jet_etagi, jet_mtagi, jet_btagi, jet_phii, jet_pti, jet_ei, jet_qi = [], [], [],[], [], [], [], [], []
    
    jet_maski, jet_etai, jet_etagi, jet_mtagi, jet_btagi, jet_phii, jet_pti, jet_ei, jet_qi = [], [], [],[], [], [], [], [], []
    nJets_ii.append(len(event.jet_eta))
        
    # 2LSS: q11, q12, q21, q22, b1, b2, b3, b4, l3, l4
    
    
    # iter_jet_phi = event.fm_phys_phi
    # iter_jet_eta = event.fm_phys_eta 
    # iter_jet_pt = event.fm_phys_pt
    # iter_jet_e = event.fm_phys_e
    # iter_jet_btag = np.zeros(len(event.fm_phys_phi))
  
    iter_jet_eta = event.jet_eta
    iter_jet_phi = event.jet_phi 
    iter_jet_pt =event.jet_pt_NOSYS
    iter_jet_e = event.jet_e_NOSYS
    iter_jet_btag = event.jet_GN2v01_FixedCutBEff_85_select
    
    maxjets = 16
    for i in range(event.nElectrons):
        lepton_etai.append(event.el_eta[i])
        lepton_etagi.append(1)
        lepton_mtagi.append(0)
        lepton_btagi.append(0)
        lepton_maski.append(True)
        lepton_phii.append(event.el_phi[i])
        lepton_pti.append(event.el_pt_NOSYS[i])
        lepton_ei.append(event.el_e_NOSYS[i])
        # lepton_ei.append(event.el_e[i])
        lepton_qi.append(event.el_charge[i])
    for i in range(event.nMuons):
        lepton_etai.append(event.mu_eta[i])
        lepton_etagi.append(0)
        lepton_mtagi.append(1)
        lepton_btagi.append(0)
        lepton_maski.append(True)
        lepton_phii.append(event.mu_phi[i])
        lepton_pti.append(event.mu_pt_NOSYS[i])
        lepton_ei.append(event.mu_e_NOSYS[i])
        # lepton_ei.append(event.mu_e[i])
        lepton_qi.append(event.mu_charge[i])
    for j in range(maxjets):
        if j >= len(iter_jet_eta):
            jet_maski.append(False)
            jet_etai.append(0)
            jet_btagi.append(0)
            jet_phii.append(0)
            jet_pti.append(0)
            jet_ei.append(0)
        else:
            jet_maski.append(True)
            jet_etai.append(iter_jet_eta[j])
            jet_btagi.append(int(iter_jet_btag[j]==chr(1)))
            jet_phii.append(iter_jet_phi[j])
            jet_pti.append(iter_jet_pt[j])
            jet_ei.append(iter_jet_e[j])
            # jet_ei.append(event.jet_e[j])
            
    met_met.append(event.met_met_NOSYS)
    met_sin_phi.append(np.cos(event.met_phi_NOSYS))
    met_cos_phi.append(np.sin(event.met_phi_NOSYS))
    # met_significance.append(event.met_significance_NOSYS)

    jet_eta.append(jet_etai)
    jet_etag.append(jet_etagi)
    jet_mtag.append(jet_mtagi)
    jet_mask.append(jet_maski)
    jet_btag.append(jet_btagi)
    jet_phi.append(jet_phii)
    jet_pt.append(jet_pti)
    jet_e.append(jet_ei)
    jet_q.append(jet_qi)

    lepton_eta.append(lepton_etai)
    lepton_etag.append(lepton_etagi)
    lepton_mtag.append(lepton_mtagi)
    lepton_mask.append(lepton_maski)
    lepton_btag.append(lepton_btagi)
    lepton_phi.append(lepton_phii)
    lepton_pt.append(lepton_pti)
    lepton_e.append(lepton_ei)
    lepton_q.append(lepton_qi)
    
    
    
    
    # After we obtain all the objects, we can perform truth matching when the sample is indeed four top. Otherwise we can just continue.
    
    if not (addTruthInformation):
        continue

    
    b_recoj_index  = event.b_recoj_index
    wd1_recoj_index = event.wd1_recoj_index
    wd2_recoj_index = event.wd2_recoj_index
    wd2_recoj_index = event.wd2_recoj_index
    wd1_el_index =event.wd1_el_index
    wd1_mu_index = event.wd1_mu_index
    wd2_el_index =event.wd2_el_index
    wd2_mu_index = event.wd2_mu_index

    all_objects_eta = np.concatenate([np.asarray(jet_etai), np.asarray(lepton_etai)])
    all_objects_phi = np.concatenate([np.asarray(jet_phii), np.asarray(lepton_phii)])    

    if do_fuzzy_matching:
        iter_objects_phi = all_objects_phi
        iter_objects_eta = all_objects_eta
    
        b_recoj_index_new = [-1,-1,-1,-1]
        Wd1_recoj_index_new = [-1,-1,-1,-1]
        Wd2_recoj_index_new = [-1,-1,-1,-1]

        b_recoj_deltaR_new = [-1,-1,-1,-1]
        Wd1_recoj_deltaR_new = [-1,-1,-1,-1]
        Wd2_recoj_deltaR_new = [-1,-1,-1,-1]
        
        wd1_el_index_new = [-1,-1,-1,-1]
        wd2_el_index_new = [-1,-1,-1,-1]
        wd1_mu_index_new = [-1,-1,-1,-1]
        wd2_mu_index_new = [-1,-1,-1,-1]

        minimum_delta_r_tolerance_b  = 0.4
        minimum_delta_r_tolerance_q = 0.4
        match_partons_to_jets(
            event.parton_b_eta,
            event.parton_b_phi,
            None,  # no pdgId filtering
            iter_objects_eta,
            iter_objects_phi,
            b_recoj_index_new,
            b_recoj_deltaR_new,
            minimum_delta_r_tolerance_b,
            None,
            None,
            None,
            event.b_recoj_index,
            [maximum_deltaR_bparton],
            label="b"
        )

        if maxjets in b_recoj_index_new:
            print(b_recoj_index_new)
            exit()
            
        # Match W decay 1
        match_partons_to_jets(
            event.parton_Wdecay1_eta,
            event.parton_Wdecay1_phi,
            event.parton_Wdecay1_pdgId,
            iter_objects_eta,
            iter_objects_phi,
            Wd1_recoj_index_new,
            Wd1_recoj_deltaR_new,
            minimum_delta_r_tolerance_q,
            wd1_el_index_new,
            wd1_mu_index_new, 
            event.nElectrons,
            event.wd1_recoj_index,
            [maximum_deltaR_Wd1],
            skip_neutrinos=True,
            label="Wd1"
        )

        # Match W decay 2
        match_partons_to_jets(
            event.parton_Wdecay2_eta,
            event.parton_Wdecay2_phi,
            event.parton_Wdecay2_pdgId,
            iter_objects_eta,
            iter_objects_phi,
            Wd2_recoj_index_new,
            Wd2_recoj_deltaR_new,
            minimum_delta_r_tolerance_q,
            wd2_el_index_new,
            wd2_mu_index_new ,
            event.nElectrons,
            event.wd2_recoj_index,

            [maximum_deltaR_Wd2],
            skip_neutrinos=True,
            label="Wd2"
        )
        
        if doDEBUG:
                print("before fixes")
                print(b_recoj_index_new)
                print(Wd1_recoj_index_new)
                print(Wd2_recoj_index_new)
            
        # Rid of the cases where there are multiple matches
        resolve_internal_duplicates(Wd1_recoj_index_new, Wd1_recoj_deltaR_new)
        resolve_internal_duplicates(Wd2_recoj_index_new, Wd2_recoj_deltaR_new)
        resolve_internal_duplicates(b_recoj_index_new, b_recoj_deltaR_new)       
        resolve_conflict(b_recoj_index_new, b_recoj_deltaR_new, Wd1_recoj_index_new, Wd1_recoj_deltaR_new)
        resolve_conflict(b_recoj_index_new, b_recoj_deltaR_new, Wd2_recoj_index_new, Wd2_recoj_deltaR_new)
        resolve_conflict(Wd1_recoj_index_new, Wd1_recoj_deltaR_new, Wd2_recoj_index_new, Wd2_recoj_deltaR_new)
        
        
        if doDEBUG:
            print("post fixes")
            print(b_recoj_index_new, " -> ", event.b_recoj_index)
            print(Wd1_recoj_index_new," -> ", event.wd1_recoj_index)
            print(Wd2_recoj_index_new," -> ", event.wd2_recoj_index)
        
        # print(event.wd1_el_index)
        # print(event.wd1_mu_index)
        # print(event.wd2_el_index)
        # print(event.wd2_mu_index)
        
        # print(wd1_el_index_new)
        # print(wd1_mu_index_new)
        # print(wd2_el_index_new)
        # print(wd2_mu_index_new)
        
        
    
        ## Time to get the leptons:
        b_count = 4 - b_recoj_index_new.count(-1)
        Wd1_count = 4-  Wd1_recoj_index_new.count(-1)
        Wd2_count = 4 - Wd2_recoj_index_new.count(-1)

        sum_count = b_count + Wd1_count+ Wd2_count
        
        running_count.append(sum_count)
        if sum_count == 11:
            exit()

    first_lt = 0
    first_ht = 0
    
    iter_b_recoj_index = event.b_recoj_index
    iter_wd1_recoj_index = event.wd1_recoj_index
    iter_wd2_recoj_index = event.wd2_recoj_index
    iter_wd1_el_index = event.wd1_el_index
    iter_wd2_el_index = event.wd2_el_index
    iter_wd1_mu_index = event.wd1_mu_index
    iter_wd2_mu_index = event.wd2_mu_index
    
    # print(max(wd1_el_index_new))
    # print(max(wd2_el_index_new))
    # print(max(wd1_mu_index_new))
    # print(max(wd2_mu_index_new))
    # if(3 in wd1_el_index_new):
    #     print("found 3, exit")
    #     exit()
    # if(3 in wd2_el_index_new):
    #     print("found 3, exit")
    #     exit()
    # if(3 in wd1_mu_index_new):
    #     print("found 3, exit")
    #     exit()
    # if(3 in wd2_mu_index_new):
    #     print("found 3, exit")
    #     exit()
    if doDeltaRMatching: 
        iter_b_recoj_index = b_recoj_index_new
        iter_wd1_recoj_index = Wd1_recoj_index_new
        iter_wd2_recoj_index = Wd2_recoj_index_new    
        iter_wd1_el_index = wd1_el_index_new
        iter_wd2_el_index = wd2_el_index_new
        iter_wd1_mu_index = wd1_mu_index_new
        iter_wd2_mu_index = wd2_mu_index_new
        
        
    # This is where be match the event to the four top geometry.
#################################### 2LSS ######################################## 
    if do2LSS:
        for t in range(4):
            if event.parton_top_isHadronic[t] == 1:
                if first_ht == 0:
                    b1.append(iter_b_recoj_index[t] if iter_b_recoj_index[t] != -1 else -1)
                    q11.append(iter_wd1_recoj_index[t] if iter_wd1_recoj_index[t] != -1 else -1)
                    q12.append(iter_wd2_recoj_index[t] if iter_wd2_recoj_index[t] != -1 else -1)

                    #check for jet matching errors
                    if iter_b_recoj_index[t] == iter_wd1_recoj_index[t] and iter_wd1_recoj_index[t] != -1:
                        b1[-1] = -1
                        q11[-1] = -1
                    if iter_wd1_recoj_index[t] == iter_wd2_recoj_index[t] and iter_wd1_recoj_index[t] != -1:
                        q11[-1] = -1
                        q12[-1] = -1
                    if iter_b_recoj_index[t] == iter_wd2_recoj_index[t] and iter_wd2_recoj_index[t] != -1:
                        b1[-1] = -1
                        q12[-1] = -1
                    first_ht = 1
                else:
                    b2.append(iter_b_recoj_index[t] if iter_b_recoj_index[t] != -1 else -1)
                    q21.append(iter_wd1_recoj_index[t] if iter_wd1_recoj_index[t] != -1 else -1)
                    q22.append(iter_wd2_recoj_index[t] if iter_wd2_recoj_index[t] != -1 else -1)
                    
                    # check for matching errors
                    if iter_b_recoj_index[t] == iter_wd1_recoj_index[t] and iter_wd1_recoj_index[t] != -1:
                        b2[-1] = -1
                        q21[-1] = -1
                    if iter_wd1_recoj_index[t] == iter_wd2_recoj_index[t] and iter_wd1_recoj_index[t] != -1:
                        q21[-1] = -1
                        q22[-1] = -1
                    if iter_b_recoj_index[t] == iter_wd2_recoj_index[t] and iter_wd2_recoj_index[t] != -1:
                        b2[-1] = -1
                        q22[-1] = -1
                        
            else:
                leptonic_W_d1_pdg_id.append(event.parton_Wdecay1_pdgId[t])
                leptonic_W_d2_pdg_id.append(event.parton_Wdecay2_pdgId[t])
                if first_lt == 0:
                    b3.append(iter_b_recoj_index[t] if iter_b_recoj_index[t] != -1 else -1)
                    ## Find the neutrino
                    if (np.abs(event.parton_Wdecay1_pdgId[t]) == 12 or np.abs(event.parton_Wdecay1_pdgId[t]) == 14): 
                        W_leptonic_nu3_truth_m.append(event.parton_Wdecay1_m[t])
                        W_leptonic_nu3_truth_pt.append(event.parton_Wdecay1_pt[t])
                        W_leptonic_nu3_truth_eta.append(event.parton_Wdecay1_eta[t])
                        W_leptonic_nu3_truth_phi.append(event.parton_Wdecay1_phi[t])
                        W_leptonic_lep3_truth_m.append(event.parton_Wdecay2_m[t])
                        W_leptonic_lep3_truth_pt.append(event.parton_Wdecay2_pt[t])
                        W_leptonic_lep3_truth_eta.append(event.parton_Wdecay2_eta[t])
                        W_leptonic_lep3_truth_phi.append(event.parton_Wdecay2_phi[t])
                    else:
                        W_leptonic_nu3_truth_m.append(event.parton_Wdecay2_m[t])
                        W_leptonic_nu3_truth_pt.append(event.parton_Wdecay2_pt[t])
                        W_leptonic_nu3_truth_eta.append(event.parton_Wdecay2_eta[t])
                        W_leptonic_nu3_truth_phi.append(event.parton_Wdecay2_phi[t])
                        W_leptonic_lep3_truth_m.append(event.parton_Wdecay1_m[t])
                        W_leptonic_lep3_truth_pt.append(event.parton_Wdecay1_pt[t])
                        W_leptonic_lep3_truth_eta.append(event.parton_Wdecay1_eta[t])
                        W_leptonic_lep3_truth_phi.append(event.parton_Wdecay1_phi[t])
                    if iter_wd1_el_index[t] != -1:
                        l3.append(iter_wd1_el_index[t])
                    elif iter_wd2_el_index[t] != -1:
                        l3.append(iter_wd2_el_index[t])
                    elif iter_wd1_mu_index[t] != -1:
                        l3.append(iter_wd1_mu_index[t] )
                    elif iter_wd2_mu_index[t] != -1:
                        l3.append(iter_wd2_mu_index[t] )
                    else:
                        l3.append(-1)
                    first_lt = 1
                else:
                    b4.append(iter_b_recoj_index[t] if iter_b_recoj_index[t] != -1 else -1)
                    if (np.abs(event.parton_Wdecay1_pdgId[t]) == 12 or np.abs(event.parton_Wdecay1_pdgId[t]) == 14): 
                        W_leptonic_nu4_truth_m.append(event.parton_Wdecay1_m[t])
                        W_leptonic_nu4_truth_pt.append(event.parton_Wdecay1_pt[t])
                        W_leptonic_nu4_truth_eta.append(event.parton_Wdecay1_eta[t])
                        W_leptonic_nu4_truth_phi.append(event.parton_Wdecay1_phi[t])
                        W_leptonic_lep4_truth_m.append(event.parton_Wdecay2_m[t])
                        W_leptonic_lep4_truth_pt.append(event.parton_Wdecay2_pt[t])
                        W_leptonic_lep4_truth_eta.append(event.parton_Wdecay2_eta[t])
                        W_leptonic_lep4_truth_phi.append(event.parton_Wdecay2_phi[t])
                    else:
                        W_leptonic_nu4_truth_m.append(event.parton_Wdecay2_m[t])
                        W_leptonic_nu4_truth_pt.append(event.parton_Wdecay2_pt[t])
                        W_leptonic_nu4_truth_eta.append(event.parton_Wdecay2_eta[t])
                        W_leptonic_nu4_truth_phi.append(event.parton_Wdecay2_phi[t])
                        W_leptonic_lep4_truth_m.append(event.parton_Wdecay1_m[t])
                        W_leptonic_lep4_truth_pt.append(event.parton_Wdecay1_pt[t])
                        W_leptonic_lep4_truth_eta.append(event.parton_Wdecay1_eta[t])
                        W_leptonic_lep4_truth_phi.append(event.parton_Wdecay1_phi[t])
                        
                    if iter_wd1_el_index[t] != -1:
                        l4.append(iter_wd1_el_index[t])
                    elif iter_wd2_el_index[t] != -1:
                        l4.append(iter_wd2_el_index[t])
                    elif iter_wd1_mu_index[t] != -1:
                        l4.append(iter_wd1_mu_index[t])# + event.nElectrons)
                    elif iter_wd2_mu_index[t] != -1:
                        l4.append(iter_wd2_mu_index[t]) #+ event.nElectrons)
                    else:
                        l4.append(-1)
#################################### 2LSS ######################################## 

#################################### 3L ######################################## 

        for t in range(4):
            if event.parton_top_isHadronic[t] == 1:
                # Hadronic Top
                b_Hadronic.append(iter_b_recoj_index[t] if iter_b_recoj_index[t] != -1 else -1)
                q1_Hadronic.append(iter_wd1_recoj_index[t] if iter_wd1_recoj_index[t] != -1 else -1)
                q2_Hadronic.append(iter_wd2_recoj_index[t] if iter_wd2_recoj_index[t] != -1 else -1)
            else:
                # Leptonic Top
                b_Leptonic.append(iter_b_recoj_index[t] if iter_b_recoj_index[t] != -1 else -1)
                l_leptonic.append()
                
                if (np.abs(event.parton_Wdecay1_pdgId[t]) == 12 or np.abs(event.parton_Wdecay1_pdgId[t]) == 14): 
                        W_leptonic_nu4_truth_m.append(event.parton_Wdecay1_m[t])
                        W_leptonic_nu4_truth_pt.append(event.parton_Wdecay1_pt[t])
                        W_leptonic_nu4_truth_eta.append(event.parton_Wdecay1_eta[t])
                        W_leptonic_nu4_truth_phi.append(event.parton_Wdecay1_phi[t])
                        W_leptonic_lep4_truth_m.append(event.parton_Wdecay2_m[t])
                        W_leptonic_lep4_truth_pt.append(event.parton_Wdecay2_pt[t])
                        W_leptonic_lep4_truth_eta.append(event.parton_Wdecay2_eta[t])
                        W_leptonic_lep4_truth_phi.append(event.parton_Wdecay2_phi[t])
                    else:    
                        W_leptonic_nu4_truth_m.append(event.parton_Wdecay2_m[t])
                        W_leptonic_nu4_truth_pt.append(event.parton_Wdecay2_pt[t])
                        W_leptonic_nu4_truth_eta.append(event.parton_Wdecay2_eta[t])
                        W_leptonic_nu4_truth_phi.append(event.parton_Wdecay2_phi[t])
                        W_leptonic_lep4_truth_m.append(event.parton_Wdecay1_m[t])
                        W_leptonic_lep4_truth_pt.append(event.parton_Wdecay1_pt[t])
                        W_leptonic_lep4_truth_eta.append(event.parton_Wdecay1_eta[t])
                        W_leptonic_lep4_truth_phi.append(event.parton_Wdecay1_phi[t])
                        
                        
                        
################################### 2LSS ######################################## 
    if do2LSS:
        for t in range(4):
            if event.parton_top_isHadronic[t] == 1:
                if first_ht == 0:
                    b1.append(iter_b_recoj_index[t] if iter_b_recoj_index[t] != -1 else -1)
                    q11.append(iter_wd1_recoj_index[t] if iter_wd1_recoj_index[t] != -1 else -1)
                    q12.append(iter_wd2_recoj_index[t] if iter_wd2_recoj_index[t] != -1 else -1)

                    #check for jet matching errors
                    if iter_b_recoj_index[t] == iter_wd1_recoj_index[t] and iter_wd1_recoj_index[t] != -1:
                        b1[-1] = -1
                        q11[-1] = -1
                    if iter_wd1_recoj_index[t] == iter_wd2_recoj_index[t] and iter_wd1_recoj_index[t] != -1:
                        q11[-1] = -1
                        q12[-1] = -1
                    if iter_b_recoj_index[t] == iter_wd2_recoj_index[t] and iter_wd2_recoj_index[t] != -1:
                        b1[-1] = -1
                        q12[-1] = -1
                    first_ht = 1
                else:
                    b2.append(iter_b_recoj_index[t] if iter_b_recoj_index[t] != -1 else -1)
                    q21.append(iter_wd1_recoj_index[t] if iter_wd1_recoj_index[t] != -1 else -1)
                    q22.append(iter_wd2_recoj_index[t] if iter_wd2_recoj_index[t] != -1 else -1)
                    
                    # check for matching errors
                    if iter_b_recoj_index[t] == iter_wd1_recoj_index[t] and iter_wd1_recoj_index[t] != -1:
                        b2[-1] = -1
                        q21[-1] = -1
                    if iter_wd1_recoj_index[t] == iter_wd2_recoj_index[t] and iter_wd1_recoj_index[t] != -1:
                        q21[-1] = -1
                        q22[-1] = -1
                    if iter_b_recoj_index[t] == iter_wd2_recoj_index[t] and iter_wd2_recoj_index[t] != -1:
                        b2[-1] = -1
                        q22[-1] = -1
                        
            else:
                leptonic_W_d1_pdg_id.append(event.parton_Wdecay1_pdgId[t])
                leptonic_W_d2_pdg_id.append(event.parton_Wdecay2_pdgId[t])
                if first_lt == 0:
                    b3.append(iter_b_recoj_index[t] if iter_b_recoj_index[t] != -1 else -1)
                    ## Find the neutrino
                    if (np.abs(event.parton_Wdecay1_pdgId[t]) == 12 or np.abs(event.parton_Wdecay1_pdgId[t]) == 14): 
                        W_leptonic_nu3_truth_m.append(event.parton_Wdecay1_m[t])
                        W_leptonic_nu3_truth_pt.append(event.parton_Wdecay1_pt[t])
                        W_leptonic_nu3_truth_eta.append(event.parton_Wdecay1_eta[t])
                        W_leptonic_nu3_truth_phi.append(event.parton_Wdecay1_phi[t])
                        W_leptonic_lep3_truth_m.append(event.parton_Wdecay2_m[t])
                        W_leptonic_lep3_truth_pt.append(event.parton_Wdecay2_pt[t])
                        W_leptonic_lep3_truth_eta.append(event.parton_Wdecay2_eta[t])
                        W_leptonic_lep3_truth_phi.append(event.parton_Wdecay2_phi[t])
                    else:
                        W_leptonic_nu3_truth_m.append(event.parton_Wdecay2_m[t])
                        W_leptonic_nu3_truth_pt.append(event.parton_Wdecay2_pt[t])
                        W_leptonic_nu3_truth_eta.append(event.parton_Wdecay2_eta[t])
                        W_leptonic_nu3_truth_phi.append(event.parton_Wdecay2_phi[t])
                        W_leptonic_lep3_truth_m.append(event.parton_Wdecay1_m[t])
                        W_leptonic_lep3_truth_pt.append(event.parton_Wdecay1_pt[t])
                        W_leptonic_lep3_truth_eta.append(event.parton_Wdecay1_eta[t])
                        W_leptonic_lep3_truth_phi.append(event.parton_Wdecay1_phi[t])
                    if iter_wd1_el_index[t] != -1:
                        l3.append(iter_wd1_el_index[t])
                    elif iter_wd2_el_index[t] != -1:
                        l3.append(iter_wd2_el_index[t])
                    elif iter_wd1_mu_index[t] != -1:
                        l3.append(iter_wd1_mu_index[t] )
                    elif iter_wd2_mu_index[t] != -1:
                        l3.append(iter_wd2_mu_index[t] )
                    else:
                        l3.append(-1)
                    first_lt = 1
                else:
                    b4.append(iter_b_recoj_index[t] if iter_b_recoj_index[t] != -1 else -1)
                    if (np.abs(event.parton_Wdecay1_pdgId[t]) == 12 or np.abs(event.parton_Wdecay1_pdgId[t]) == 14): 
                        W_leptonic_nu4_truth_m.append(event.parton_Wdecay1_m[t])
                        W_leptonic_nu4_truth_pt.append(event.parton_Wdecay1_pt[t])
                        W_leptonic_nu4_truth_eta.append(event.parton_Wdecay1_eta[t])
                        W_leptonic_nu4_truth_phi.append(event.parton_Wdecay1_phi[t])
                        W_leptonic_lep4_truth_m.append(event.parton_Wdecay2_m[t])
                        W_leptonic_lep4_truth_pt.append(event.parton_Wdecay2_pt[t])
                        W_leptonic_lep4_truth_eta.append(event.parton_Wdecay2_eta[t])
                        W_leptonic_lep4_truth_phi.append(event.parton_Wdecay2_phi[t])
                    else:    
                        W_leptonic_nu4_truth_m.append(event.parton_Wdecay2_m[t])
                        W_leptonic_nu4_truth_pt.append(event.parton_Wdecay2_pt[t])
                        W_leptonic_nu4_truth_eta.append(event.parton_Wdecay2_eta[t])
                        W_leptonic_nu4_truth_phi.append(event.parton_Wdecay2_phi[t])
                        W_leptonic_lep4_truth_m.append(event.parton_Wdecay1_m[t])
                        W_leptonic_lep4_truth_pt.append(event.parton_Wdecay1_pt[t])
                        W_leptonic_lep4_truth_eta.append(event.parton_Wdecay1_eta[t])
                        W_leptonic_lep4_truth_phi.append(event.parton_Wdecay1_phi[t])
                        
                    if iter_wd1_el_index[t] != -1:
                        l4.append(iter_wd1_el_index[t])
                    elif iter_wd2_el_index[t] != -1:
                        l4.append(iter_wd2_el_index[t])
                    elif iter_wd1_mu_index[t] != -1:
                        l4.append(iter_wd1_mu_index[t])# + event.nElectrons)
                    elif iter_wd2_mu_index[t] != -1:
                        l4.append(iter_wd2_mu_index[t]) #+ event.nElectrons)
                    else:
                        l4.append(-1)
#################################### 3L ######################################## 


hf = h5py.File(name, 'w')
inputs = hf.create_group('INPUTS')

jet = inputs.create_group('Jet')
jet_mask_data = jet.create_dataset('MASK', data=jet_mask)
jet_eta_data = jet.create_dataset('eta', data=jet_eta)
jet_btag_data = jet.create_dataset('btag', data=jet_btag)
jet_phi_data = jet.create_dataset('phi', data=jet_phi)
jet_pt_data = jet.create_dataset('pt', data=jet_pt)
jet_e_data = jet.create_dataset('e', data=jet_e)

lepton = inputs.create_group('Lepton')
lepton_mask_data = lepton.create_dataset('MASK', data=lepton_mask)
lepton_eta_data = lepton.create_dataset('eta', data=lepton_eta)
lepton_btag_data = lepton.create_dataset('btag', data=lepton_btag)
lepton_phi_data = lepton.create_dataset('phi', data=lepton_phi)
lepton_pt_data = lepton.create_dataset('pt', data=lepton_pt)
lepton_e_data = lepton.create_dataset('e', data=lepton_e)
lepton_q_data = lepton.create_dataset('q', data=lepton_q)
lepton_etag_data = lepton.create_dataset('etag', data=lepton_etag)
lepton_mtag_data = lepton.create_dataset('mtag', data=lepton_mtag)


if not (addTruthInformation):
    exit()


print(maximum_deltaR_bparton)  
print(maximum_deltaR_Wd1)  
print(maximum_deltaR_Wd2)  

print("q11: " , min(q11) , max(q11))
print("q12: " ,min(q12) , max(q12))
print("q21: ",min(q21) , max(q21))
print("q22: ",min(q22) , max(q22))
print("l3: ",min(l3) , max(l3))
print("l4: ",min(l4) , max(l4))
print("b1: ",min(b1) , max(b1))
print("b2: ",min(b2) , max(b2))
print("b3: ",min(b3) , max(b3))
print("b4: ",min(b4) , max(b4))
# for i in range(30):
#     print(leptonic_W_d1_pdg_id[i],leptonic_W_d2_pdg_id[i],W_leptonic_nu4_truth_m[i],W_leptonic_lep4_truth_m[i])
    
    
plotDirectory = "matching_plots"
os.makedirs(plotDirectory, exist_ok=True)

## Plots all the matching for each of the objects. If it is -1, it means there were not matching
plt.figure()
plt.hist(b1, bins=50, range=(-2, 10),  histtype='step', color='blue', label="b1") 
plt.savefig(plotDirectory+"/b1.pdf") #Save the plot

plt.figure()
plt.hist(b2, bins=50, range=(-2, 10),  histtype='step', color='blue', label="b1") 
plt.savefig(plotDirectory + "/b2.pdf") #Save the plot

plt.figure()
plt.hist(b3, bins=50, range=(-2, 10),  histtype='step', color='blue', label="b1") 
plt.savefig(plotDirectory +"/b3.pdf") #Save the plot

plt.figure()
plt.hist(b4, bins=50, range=(-2, 10),  histtype='step', color='blue', label="b1")
plt.savefig(plotDirectory +"/b4.pdf") #Save the plot

plt.figure()
plt.hist(q11, bins=50, range=(-2, 10),  histtype='step', color='blue', label="b1")
plt.savefig(plotDirectory +"/q11.pdf") #Save the plot

plt.figure()
plt.hist(q12, bins=50, range=(-2, 10),  histtype='step', color='blue', label="b1") 
plt.savefig(plotDirectory +"/q12.pdf") #Save the plot
#####################################################################################################################

assigned_b_jet_vect = []
assigned_w_had_vect = []
assigned_w_lep_vect = []
assigned_w_total_vect = []
assigned_t_had_vect = []
assigned_t_lep_vect = []
assigned_t_total_vect = []
assigned_q_jet_vect = []
assigned_jet_vect = []
assigned_objects = []

print(len(b1))
for ib1, ib2,ib3,ib4,iq11,iq12,iq21,iq22,il3,il4 in zip(b1, b2,b3,b4,q11,q12,q21,q22,l3,l4):
    assigned_b_jet = 0
    assigned_w_had = 0
    assigned_w_lep = 0
    assigned_w_total = 0
    assigned_t_had = 0
    assigned_t_lep = 0
    assigned_t_total = 0
    assigned_q_jet= 0 
    assigned_lep = 0

    if (ib1 !=-1):
        assigned_b_jet = assigned_b_jet  +1
    if (ib2 !=-1):
        assigned_b_jet = assigned_b_jet  +1
    if (ib3 !=-1):
        assigned_b_jet = assigned_b_jet  +1
    if (ib4 !=-1):
        assigned_b_jet = assigned_b_jet  +1
    if (iq11 !=-1):
        assigned_q_jet = assigned_q_jet  +1
    if (iq12 !=-1):
        assigned_q_jet = assigned_q_jet  +1
    if (iq21 !=-1):
        assigned_q_jet = assigned_q_jet  +1
    if (iq22 !=-1):
        assigned_q_jet = assigned_q_jet  +1
    if (il3 !=-1):
        assigned_lep = assigned_lep  +1
    if (il4 !=-1):
        assigned_lep = assigned_lep  +1    


    if((iq11 != -1) and (iq12 != -1)):
        assigned_w_had = assigned_w_had +1
        if (ib1 !=-1):
            assigned_t_had = assigned_t_had +1
    if((iq21 != -1) and (iq22 != -1)):
        assigned_w_had = assigned_w_had +1  
        if (ib2 !=-1):
            assigned_t_had = assigned_t_had +1             
    if(il3 != -1):
        assigned_w_lep = assigned_w_lep +1
        if (ib3 !=-1):
            assigned_t_lep = assigned_t_lep  +1      
    if(il4 != -1):
        assigned_w_lep = assigned_w_lep +1    
        if (ib4 !=-1):
            assigned_t_lep = assigned_t_lep  +1          

    assigned_w_had_vect.append(assigned_w_had)
    assigned_w_lep_vect.append(assigned_w_lep)
    assigned_w_total_vect.append(assigned_w_had +assigned_w_lep)
    assigned_t_had_vect.append(assigned_t_had)
    assigned_t_lep_vect.append(assigned_t_lep)
    assigned_t_total_vect.append(assigned_t_had +assigned_t_lep)                     
    assigned_b_jet_vect.append(assigned_b_jet)
    assigned_q_jet_vect.append(assigned_q_jet)
    assigned_jet_vect.append(assigned_q_jet + assigned_b_jet)
    assigned_objects.append(assigned_q_jet + assigned_b_jet+assigned_lep)

print(len(assigned_b_jet_vect))
plt.figure()    
plt.hist(assigned_b_jet_vect, bins=12, range=(-1, 5),  histtype='step', color='blue', label="b1") #histtype='step' for a line plot
plt.yscale('log')
plt.ylim(1, len(assigned_t_had_vect))  # Set y-axis range from 0 to 100,000
plt.savefig(plotDirectory +"/correctly_assigned_bjets.pdf") #Save the plot

plt.figure()    
plt.hist(assigned_w_had_vect, bins=12, range=(-1, 5),  histtype='step', color='blue', label="b1") #histtype='step' for a line plot
plt.savefig(plotDirectory +"/correctly_assigned_w_had.pdf") #Save the plot

plt.figure()    
plt.hist(assigned_w_lep_vect, bins=12, range=(-1, 5),  histtype='step', color='blue', label="b1") #histtype='step' for a line plot
plt.savefig(plotDirectory +"/correctly_assigned_w_lep.pdf") #Save the plot

plt.figure()    
plt.hist(assigned_w_total_vect, bins=12, range=(-1, 5),  histtype='step', color='blue', label="b1") #histtype='step' for a line plot
plt.savefig(plotDirectory +"/correctly_assigned_w_total.pdf") #Save the plot

plt.figure()    
plt.hist(assigned_t_had_vect, bins=12, range=(-1, 5),  histtype='step', color='blue', label="b1") #histtype='step' for a line plot
plt.yscale('log')
plt.ylim(1, len(assigned_t_had_vect))  # Set y-axis range from 0 to 100,000
plt.savefig(plotDirectory +"/correctly_assigned_t_had.pdf") #Save the plot

plt.figure()    
plt.hist(assigned_t_lep_vect, bins=12, range=(-1, 5),  histtype='step', color='blue', label="b1") #histtype='step' for a line plot
plt.yscale('log')
plt.ylim(1, len(assigned_t_had_vect))  # Set y-axis range from 0 to 100,000
plt.savefig(plotDirectory +"/correctly_assigned_t_lep.pdf") #Save the plot

plt.figure()    
plt.hist(assigned_t_total_vect, bins=12, range=(-1, 5),  histtype='step', color='blue', label="b1") #histtype='step' for a line plot
plt.yscale('log')
plt.ylim(1, len(assigned_b_jet_vect))  # Set y-axis range from 0 to 100,000
plt.savefig(plotDirectory +"/correctly_assigned_t_total.pdf") #Save the plot

plt.figure()    
plt.hist(assigned_q_jet_vect, bins=12, range=(-1, 5),  histtype='step', color='blue', label="b1") #histtype='step' for a line plot
plt.yscale('log')
plt.ylim(1, len(assigned_t_had_vect))  # Set y-axis range from 0 to 100,000
plt.savefig(plotDirectory +"/correctly_assigned_qjets.pdf") #Save the plot

plt.figure()    
n, bins, patches = plt.hist(assigned_jet_vect, bins=9, range=(-1, 8),  histtype='step', color='blue', label="b1") #histtype='step' for a line plot
plt.yscale('log')
plt.ylim(1, len(assigned_b_jet_vect))  # Set y-axis range from 0 to 100,000
print("Bin counts:", n)
print("Bin edges:", bins)
plt.savefig(plotDirectory +"/correctly_assigned_jets.pdf") #Save the plot

cumsum = [0]  # Initialize with 0
cumsum = np.cumsum(n)
print("Cumulative sum:", cumsum)
cumsum /= cumsum[-1]  # Normalize by the last value

cumsum_reversed = cumsum
bins = [0,1,2,3,4,5,6,7,8,9]
print("Normalized Bin counts:", n)
print("Bin edges:", bins)
print("Normalized Cumulative sum:", cumsum)
print("Reversed Cumulative sum:", cumsum_reversed)


# # print(cumsum2)
# for count in n:
#     cumsum.append(cumsum[-1] + count)
    
plt.figure()    
plt.step(bins[:-1], cumsum_reversed, where='post', color='red', label="Cumulative Sum")
plt.axhline(y=0.5, color='red', linestyle='--', label='y=0.5')
plt.ylim(0, 1)  # Set y-axis range from 0 to 100,000
for i in range(len(bins) - 1):
    plt.text(bins[i], cumsum_reversed[i] + 0.02, f"{cumsum_reversed[i]:.2f}", 
             ha='center', fontsize=10, color='black')
plt.savefig(plotDirectory +"/correctly_assigned_jets_cumsum.pdf") #Save the plot


plt.figure()    
plt.hist(assigned_objects, bins=24, range=(-1, 11),  histtype='step', color='blue', label="b1") #histtype='step' for a line plot
plt.yscale('log')
plt.ylim(1, len(assigned_b_jet_vect))  # Set y-axis range from 0 to 100,000
plt.savefig(plotDirectory +"/correctly_assigned_objects.pdf") #Save the plot



############### FIND THE REGRESSIONS ###################
    
top_masses = []
top_njets = []
t1_tMass = []
t2_tMass = []
t1_wMass = []
t2_wMass = []
t1_njet = []
t2_njet = []
w1_njet = []
w2_njet =[]
t3_mbl = []
t4_mbl = []
t3_njet = []
t4_njet = []
bins = np.linspace(0, 300000, 100)  # 100 bins between 0 and 300000

for i in range(len(q11)):
    
    t1_quarks = []
    t2_quarks = []
    t3_quarks = []
    t4_quarks = []
    w1_quarks = []
    w2_quarks = []
    
    quarks = []
    if (q11[i] != -1 and  q11[i] <maxjets):
        quark = (jet_pt[i][q11[i]], jet_eta[i][q11[i]], jet_phi[i][q11[i]], jet_e[i][q11[i]])  # (pt, eta, phi, E)
        w1_quarks.append(quark)
        t1_quarks.append(quark)
        
    if(q12[i] != -1 and  q12[i] <maxjets):
        quark = (jet_pt[i][q12[i]], jet_eta[i][q12[i]], jet_phi[i][q12[i]], jet_e[i][q12[i]])  # (pt, eta, phi, E)
        w1_quarks.append(quark)
        t1_quarks.append(quark)     
           
    if(b1[i] != -1 and  b1[i] <maxjets):
        quark = (jet_pt[i][b1[i]], jet_eta[i][b1[i]], jet_phi[i][b1[i]], jet_e[i][b1[i]])  # (pt, eta, phi, E)
        t1_quarks.append(quark)         

    if (q21[i] != -1 and  q21[i] <maxjets):
        quark = (jet_pt[i][q21[i]], jet_eta[i][q21[i]], jet_phi[i][q21[i]], jet_e[i][q21[i]])  # (pt, eta, phi, E)
        w2_quarks.append(quark)
        t2_quarks.append(quark)
        
    if(q22[i] != -1 and  q22[i] <maxjets):
        quark = (jet_pt[i][q22[i]], jet_eta[i][q22[i]], jet_phi[i][q22[i]], jet_e[i][q22[i]])  # (pt, eta, phi, E)
        w2_quarks.append(quark)
        t2_quarks.append(quark)     
           
    if(b2[i] != -1 and  b2[i] <maxjets):
        quark = (jet_pt[i][b2[i]], jet_eta[i][b2[i]], jet_phi[i][b2[i]], jet_e[i][b2[i]])  # (pt, eta, phi, E)
        t2_quarks.append(quark)   
        
    if(b3[i] != -1 and  b3[i] <maxjets):
        quark = (jet_pt[i][b3[i]], jet_eta[i][b3[i]], jet_phi[i][b3[i]], jet_e[i][b3[i]])  # (pt, eta, phi, E)
        t3_quarks.append(quark)        
    if(l3[i] != -1 ):
        lepton = (lepton_pt[i][l3[i]], lepton_eta[i][l3[i]], lepton_phi[i][l3[i]], lepton_e[i][l3[i]])  # (pt, eta, phi, E)
        t3_quarks.append(lepton)  
              
    if(b4[i] != -1 and  b4[i] <maxjets):
        quark = (jet_pt[i][b4[i]], jet_eta[i][b4[i]], jet_phi[i][b4[i]], jet_e[i][b4[i]])  # (pt, eta, phi, E)
        t4_quarks.append(quark)        
    if(l4[i] != -1 ):
        lepton = (lepton_pt[i][l4[i]], lepton_eta[i][l4[i]], lepton_phi[i][l4[i]], lepton_e[i][l4[i]])  # (pt, eta, phi, E)
        t4_quarks.append(lepton)       
    
    
    
    # t1 mass:
    if len(t1_quarks) >= 1:
        t1_mass = calculate_invariant_mass(*t1_quarks)
        t1_tMass.append(t1_mass)
        t1_njet.append(len(t1_quarks))
    else:
        t1_tMass.append(-1)
        t1_njet.append(0)


    # w1 mass:
    if len(w1_quarks) >= 1:
        w1_mass = calculate_invariant_mass(*w1_quarks)
        t1_wMass.append(w1_mass)
        w1_njet.append(len(w1_quarks))
    else:
        t1_wMass.append(-1)
        w1_njet.append(0)
        
    # t2 mass
    if len(t2_quarks) >= 1:
        t2_mass = calculate_invariant_mass(*t2_quarks)
        t2_tMass.append(t2_mass)
        t2_njet.append(len(t2_quarks))
    else:
        t2_tMass.append(-1)
        t2_njet.append(0)
        
    if len(w2_quarks) >= 1:
        w2_mass = calculate_invariant_mass(*w2_quarks)
        t2_wMass.append(w2_mass)
        w2_njet.append(len(w2_quarks))
    else:
        t2_wMass.append(-1)
        w2_njet.append(0)
        
    # t3 mblmass
    if len(t3_quarks) == 2:
        m_bl = calculate_invariant_mass(*t3_quarks)            
        t3_mbl.append(m_bl)
        t3_njet.append(len(t3_quarks))
    else:
        # t3_mbl.append(np.nan)
        t3_mbl.append(-1)
        t3_njet.append(0)
        
    if len(t4_quarks) == 2:
        m_bl = calculate_invariant_mass(*t4_quarks)
        t4_mbl.append(m_bl)
        t4_njet.append(len(t4_quarks))
    else:
        t4_mbl.append(0)
        t4_njet.append(0)
        
# print(t3_mbl[0:20])    
# print(t3_njet[0:20])  
# exit()   
# Hadronic tops: 


plt.figure(figsize=(8,5))
t1_tMass = np.array(t1_tMass)
t2_tMass = np.array(t2_tMass)
t1_njet = np.array(t1_njet)
t2_njet = np.array(t2_njet)
w1_njet = np.array(w1_njet)
w2_njet = np.array(w2_njet)
t4_njet = np.array(t4_njet)
t3_njet = np.array(t3_njet)
t3_mbl = np.array(t3_mbl)
t4_mbl = np.array(t4_mbl)

top_masses = np.concatenate([t1_tMass, t2_tMass])
top_njets = np.concatenate([t1_njet, t2_njet])


top_masses = np.array(top_masses)
top_njets = np.array(top_njets)

plt.hist(top_masses, bins=bins, color='mediumseagreen', edgecolor='black', alpha=0.7)
# Group by njets
masses_1j = top_masses[top_njets == 1]
masses_2j = top_masses[top_njets == 2]
masses_3j = top_masses[top_njets == 3]

plt.hist([masses_1j, masses_2j, masses_3j],
         bins=bins,
         stacked=True,
         color=['red', 'blue', 'purple'],
         label=['1 jet', '2 jets', '3 jets'],
         edgecolor='black',
         alpha=0.7)


plt.axvline(x=173000, color='black', linestyle='--', linewidth=2, label='Top Mass (~173 GeV)')
plt.axvline(x=80377, color='black', linestyle='dashdot', linewidth=2, label='W Mass (~80 GeV)')
plt.title("Top Quark Mass Distribution (simulated)", fontsize=14)
plt.xlabel("Mass [GeV]", fontsize=12)
plt.ylabel("Counts", fontsize=12)
plt.grid(True, alpha=0.3)
plt.xlim(0, 300000)
plt.ylim(0, len(top_masses)/10)
plt.tight_layout()
plt.legend()
plt.savefig(plotDirectory +"/"+prefix +"_top_masses.pdf") #Save the plot

print(len(top_masses))


plt.figure(figsize=(8,5))
plt.hist(top_masses, bins=bins, color='mediumseagreen', edgecolor='black', alpha=0.7)

plt.hist([ masses_3j],
         bins=bins,
         stacked=True,
         color=[ 'purple'],
         label=['3 jets'],
         edgecolor='black',
         alpha=0.7)


plt.axvline(x=172760, color='black', linestyle='--', linewidth=2, label='Top Mass (~173 GeV)')
plt.axvline(x=80377, color='black', linestyle='dashdot', linewidth=2, label='W Mass (~80 GeV)')

plt.title("Top Quark Mass Distribution (simulated)", fontsize=14)
plt.xlabel("Mass [GeV]", fontsize=12)
plt.ylabel("Counts", fontsize=12)
plt.grid(True, alpha=0.3)
plt.xlim(0, 300000)
plt.ylim(0, len(top_masses)/10)

plt.tight_layout()
plt.legend()
plt.savefig(plotDirectory +"/"+prefix +"_topMass_just_njets3.pdf") #Save the plot


# desired_length = 2 * len(q11)
# padding_needed = desired_length - len(top_njets)
# if padding_needed > 0:
#     top_njets = np.pad(top_njets, (0, padding_needed), constant_values=0)


plot_pie_chart_from_array(top_njets,plotDirectory,prefix)


t3_nu_px,t3_nu_py, t3_nu_pz =  pt_eta_phi_to_pxyz(W_leptonic_nu3_truth_phi, W_leptonic_nu3_truth_eta, W_leptonic_nu3_truth_phi)
t4_nu_px,t4_nu_py, t4_nu_pz =  pt_eta_phi_to_pxyz(W_leptonic_nu3_truth_phi, W_leptonic_nu3_truth_eta, W_leptonic_nu3_truth_phi)


t3_nu_eta = W_leptonic_nu3_truth_eta
t3_nu_cos_phi = np.cos(W_leptonic_nu3_truth_phi)
t3_nu_sin_phi = np.sin(W_leptonic_nu3_truth_phi)

t4_nu_eta = W_leptonic_nu4_truth_eta
t4_nu_cos_phi = np.cos(W_leptonic_nu4_truth_phi)
t4_nu_sin_phi =  np.sin(W_leptonic_nu4_truth_phi)

############# SAVE EVERYTHING TO SPANNET ######################################



# Set up h5 file strutcture 
# hf = h5py.File(name, 'w')

# inputs = hf.create_group('INPUTS')
targets = hf.create_group('TARGETS')
truth = hf.create_group('TRUTH')
regressions = hf.create_group('REGRESSIONS')
# particle = hf.create_group('PARTICLE')


t1 = targets.create_group('t1')
t2 = targets.create_group('t2')
t3 = targets.create_group('t3')
t4 = targets.create_group('t4')
t1_truth = truth.create_group('t1')
t2_truth = truth.create_group('t2')
t3_truth = truth.create_group('t3')
t4_truth = truth.create_group('t4')
event_truth= truth.create_group('EVENT')

t1_regression = regressions.create_group('t1')
t2_regression = regressions.create_group('t2')
t3_regression = regressions.create_group('t3')
t4_regression = regressions.create_group('t4')

# event_regression = hf.create_group('EVENT')
event_regression = regressions.create_group('EVENT')

b1_data = t1.create_dataset('b', data=b1)
q11_data = t1.create_dataset('q1', data=q11)
q12_data = t1.create_dataset('q2', data=q12)
b2_data = t2.create_dataset('b', data=b2)
q21_data = t2.create_dataset('q1', data=q21)
q22_data = t2.create_dataset('q2', data=q22)
b3_data = t3.create_dataset('b', data=b3)
l3_data = t3.create_dataset('l', data=l3)
b4_data = t4.create_dataset('b', data=b4)
l4_data = t4.create_dataset('l', data=l4)


b1_data_truth = t1_truth.create_dataset('b', data=b1)
q11_data_truth = t1_truth.create_dataset('q1', data=q11)
q12_data_truth = t1_truth.create_dataset('q2', data=q12)
b2_data_truth = t2_truth.create_dataset('b', data=b2)
q21_data_truth = t2_truth.create_dataset('q1', data=q21)
q22_data_truth = t2_truth.create_dataset('q2', data=q22)
b3_data_truth = t3_truth.create_dataset('b', data=b3)
l3_data_truth = t3_truth.create_dataset('l', data=l3)
b4_data_truth = t4_truth.create_dataset('b', data=b4)
l4_data_truth = t4_truth.create_dataset('l', data=l4)

# event = regressions.create_group('EVENT')
# t3_nu3_eta_data = event.create_dataset('t3_neutrino_eta', data=t3_nu_eta)
# t3_nu_cos_phi_data = event.create_dataset('t3_neutrino_cos_phi', data=t3_nu_cos_phi)
# t3_nu_sin_phi_data = event.create_dataset('t3_neutrino_sin_phi', data=t3_nu_sin_phi)

# t4_nu_eta_data = event.create_dataset('t4_neutrino_eta', data=t4_nu_eta)
# t4_nu_cos_phi_data = event.create_dataset('t4_neutrino_cos_phi', data=t4_nu_cos_phi)
# t4_nu_sin_phi_data = event.create_dataset('t4_neutrino_sin_phi', data=t4_nu_sin_phi)

met = inputs.create_group('Met')
met_met_data = met.create_dataset('met_met', data=met_met)
met_cos_phi_data = met.create_dataset('met_cos_phi', data=met_cos_phi)
met_sin_phi_data = met.create_dataset('met_sin_phi', data=met_sin_phi)

# ## REGRESSIONS ##

t1_PARTICLE = t1_regression.create_group('PARTICLE')
t2_PARTICLE = t2_regression.create_group('PARTICLE')


t1_tMass_data = t1_PARTICLE.create_dataset('tMass', data=t1_tMass)
t2_tMass_data = t2_PARTICLE.create_dataset('tMass', data=t2_tMass)
t1_wMass_data = t1_regression.create_dataset('wMass', data=t1_wMass)
t2_wMass_data = t2_regression.create_dataset('wMass', data=t2_wMass)

t3_mbl_data = t3_regression.create_dataset('mbl', data=t3_mbl)
t3_nu_eta_data = t3_regression.create_dataset('nu_eta', data=t3_nu_eta)
t3_nu_cos_phi_data = t3_regression.create_dataset('nu_cos_phi', data=t3_nu_cos_phi)
t3_nu_sin_phi_data = t3_regression.create_dataset('nu_sin_phi', data=t3_nu_sin_phi)

t4_mbl_data = t4_regression.create_dataset('mbl', data=t4_mbl)
t4_nu_eta_data = t4_regression.create_dataset('nu_eta', data=t4_nu_eta)
t4_nu_cos_phi_data = t4_regression.create_dataset('nu_cos_phi', data=t4_nu_cos_phi)
t4_nu_sin_phi_data = t4_regression.create_dataset('nu_sin_phi', data=t4_nu_sin_phi)


t3_nu_eta_data = event_regression.create_dataset('t3_nu_eta', data=t3_nu_eta)
t3_nu_cos_phi_data = event_regression.create_dataset('t3_nu_cos_phi', data=t3_nu_cos_phi)
t3_nu_sin_phi_data = event_regression.create_dataset('t3_nu_sin_phi', data=t3_nu_sin_phi)

t4_nu_eta_data = event_regression.create_dataset('t4_nu_eta', data=t4_nu_eta)
t4_nu_cos_phi_data = event_regression.create_dataset('t4_nu_cos_phi', data=t4_nu_cos_phi)
t4_nu_sin_phi_data = event_regression.create_dataset('t4_nu_sin_phi', data=t4_nu_sin_phi)

t3_nu_eta_data = event_truth.create_dataset('t3_nu_eta', data=t3_nu_eta)
t3_nu_cos_phi_data = event_truth.create_dataset('t3_nu_cos_phi', data=t3_nu_cos_phi)
t3_nu_sin_phi_data = event_truth.create_dataset('t3_nu_sin_phi', data=t3_nu_sin_phi)

t4_nu_eta_data = event_truth.create_dataset('t4_nu_eta', data=t4_nu_eta)
t4_nu_cos_phi_data = event_truth.create_dataset('t4_nu_cos_phi', data=t4_nu_cos_phi)
t4_nu_sin_phi_data = event_truth.create_dataset('t4_nu_sin_phi', data=t4_nu_sin_phi)



t3_nu_px_data = event_regression.create_dataset('t3_nu_px', data=t3_nu_px)
t3_nu_py_data = event_regression.create_dataset('t3_nu_py', data=t3_nu_py)
t3_nu_pz_data = event_regression.create_dataset('t3_nu_pz', data=t3_nu_pz)
t3_nu_px_data = event_truth.create_dataset('t3_nu_px', data=t3_nu_px)
t3_nu_py_data = event_truth.create_dataset('t3_nu_py', data=t3_nu_py)
t3_nu_pz_data = event_truth.create_dataset('t3_nu_pz', data=t3_nu_pz)

t4_nu_px_data = event_regression.create_dataset('t4_nu_px', data=t4_nu_px)
t4_nu_py_data = event_regression.create_dataset('t4_nu_py', data=t4_nu_py)
t4_nu_pz_data = event_regression.create_dataset('t4_nu_pz', data=t4_nu_pz)
t4_nu_px_data = event_truth.create_dataset('t4_nu_px', data=t4_nu_px)
t4_nu_py_data = event_truth.create_dataset('t4_nu_py', data=t4_nu_py)
t4_nu_pz_data = event_truth.create_dataset('t4_nu_pz', data=t4_nu_pz)

# t2_tMass 
t2_validTops  = t2_njet ==3
t1_validTops  = t1_njet ==3
t2_validWs  = w2_njet ==2
t1_validWs  = w1_njet ==2
t3_valid_mbl  = t3_njet ==2
t4_valid_mbl  = t4_njet ==2


# exit()
for i in range(len(t2_validTops)):
    if (t2_validTops[i] == False):
        t2_tMass[i] = np.nan
    if (t1_validTops[i] == False):
        t1_tMass[i] = np.nan    
    if (t1_validWs[i] == False):
        t1_wMass[i] = np.nan 
    if (t2_validWs[i] == False):
        t2_wMass[i] = np.nan 
        
    if (t3_valid_mbl[i] == False):
        t3_mbl[i] = np.nan 
    if (t4_valid_mbl[i] == False):
        t4_mbl[i] = np.nan 
        
# print(t2_validTops)
# print(t2_tMass)
# exit()

t1_tMass_data = event_regression.create_dataset('t1_tMass', data=t1_tMass)
t2_tMass_data = event_regression.create_dataset('t2_tMass', data=t2_tMass)
t1_tMass_data = event_truth.create_dataset('t1_tMass', data=t1_tMass)
t2_tMass_data = event_truth.create_dataset('t2_tMass', data=t2_tMass)

t1_tMass_data = event_regression.create_dataset('t1_wMass', data=t1_wMass)
t2_tMass_data = event_regression.create_dataset('t2_wMass', data=t2_wMass)
t1_tMass_data = event_truth.create_dataset('t1_wMass', data=t1_wMass)
t2_tMass_data = event_truth.create_dataset('t2_wMass', data=t2_wMass)

print(t3_mbl[1:50])
print(t3_njet[1:50])

t3_mbl_data = event_regression.create_dataset('t3_mbl', data=t3_mbl)
t4_mbl_data = event_regression.create_dataset('t4_mbl', data=t4_mbl)
t3_mbl_data = event_truth.create_dataset('t3_mbl', data=t3_mbl)
t4_mbl_data = event_truth.create_dataset('t4_mbl', data=t4_mbl)


## TRUTH INFO ##
Event = inputs.create_group('Event')
Event_assigned_objects_data = Event.create_dataset('assigned_objects', data=assigned_objects)
Event_assigned_jets_data = Event.create_dataset('assigned_jets', data=assigned_jet_vect)
Event_assigned_bjets_data = Event.create_dataset('assigned_bjets', data=assigned_b_jet_vect)
Event_assigned_tops_data = Event.create_dataset('assigned_tops', data=assigned_t_total_vect)
Event_assigned_q_jets_data = Event.create_dataset('assigned_q_jets', data=assigned_q_jet_vect)
Event_assigned_had_tops_data = Event.create_dataset('assigned_had_tops', data=assigned_t_had_vect)
Event_njets_data = Event.create_dataset('nJets', data=nJets_ii)


### Do the train/test splitting uncomment this line ###
# filter_and_save_h5_random(sys.argv[1],prefix)



# with h5py.File(name, 'r') as f:
#     data = f['TARGETS']['t3']['l'][:]  # Load the entire dataset into memory
#     print(max(data))
#     print(min(data))
#     print(data)       # Print the contents    


hf.close()
