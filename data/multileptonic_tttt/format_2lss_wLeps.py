import h5py
import sys
import ROOT
import numpy as np
import glob
maxjets = 16

# Get root files from command line input

#TRAINING NTUPLES
# /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20a/*412043*/*000003* /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20d/*412043*/* /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20e/*412043*/* /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20*/*412044*/* /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20*/*700355*/* 



#412043 mc20a /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20a/*412043*/*000003*
#412043 mc20d /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20d/*412043*/* 
#412043 mc20e /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20e/*412043*/* 
#412044 mc20ade /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20*/*412044*/* 
#700355 mc20ade /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20*/*700355*/* 


# TESTING NTUPLE
#412043 mc20a /eos/atlas/atlascerngroupdisk/phys-top/topplusx/4tops2021/Run3/FullValidation/mc20a/*412043*/*000002*



chain = ROOT.TChain("reco")
folder_path = sys.argv[2]  # Make sure this is a directory, not a file
print(folder_path)
root_files = glob.glob(f"{folder_path}/*.root")
print(root_files)
# Add each file to the chain
for root_file in root_files:
    chain.Add(root_file)

print(f"Added {len(root_files)} files to the chain.")


tree = chain.CopyTree("pass_SSee_passECIDS_NOSYS||pass_SSem_passECIDS_NOSYS||pass_SSmm_NOSYS")

name = str(sys.argv[1])

if '.h5' not in name:
    name += '.h5'

# Set up h5 file strutcture 


hf = h5py.File(name, 'w')

inputs = hf.create_group('INPUTS')
targets = hf.create_group('TARGETS')
# classes = hf.create_group('CLASSIFICATIONS')

t1 = targets.create_group('t1')
t2 = targets.create_group('t2')
t3 = targets.create_group('t3')
t4 = targets.create_group('t4')


emptylist = [[0]*maxjets]*tree.GetEntries()

# jet_info
mask_jet = []
eta_jet = []
btag_jet = []
phi_jet =[]
pt_jet =[]
e_jet = []

mask_lepton = []
eta_lepton = []
etag_lepton = []
mtag_lepton = []
q_lepton = []
phi_lepton =[]
pt_lepton =[]
e_lepton = []


b1 = []
b2 = []
b3 =[]
b4 = []

q11 = []
q12 = []
q21 = []
q22 = []
l3 =[]
l4 = []

for event in tree:
    if (len(event.parton_top_isHadronic) != 4):
        continue
    if sum(event.parton_top_isHadronic) != 2:
        continue
            
        mask_jet_iter = []
        eta_jet_iter = []
        btag_jet_iter = []
        phi_jet_iter =[]
        pt_jet_iter =[]
        e_jet_iter = []

        mask_lepton_iter = []
        eta_lepton_iter = []
        etag_lepton_iter = []
        mtag_lepton_iter = []
        q_lepton_iter = []
        phi_lepton_iter =[]
        pt_lepton_iter =[]
        e_lepton_iter = []

        # first do our leptons
        for i in range(event.nElectrons):
            eta_lepton_iter.append(event.el_eta[i])
            etag_lepton_iter.append(1)
            mtag_lepton_iter.append(0)
            mask_lepton_iter.append(True)
            phi_lepton_iter.append(event.el_phi[i])
            pt_lepton_iter.append(event.el_pt_NOSYS[i])
            e_lepton_iter.append(event.el_e[i])
            q_lepton_iter.append(event.el_charge[i])
            
        for i in range(event.nMuons):
            eta_lepton_iter.append(event.mu_eta[i])
            etag_lepton_iter.append(0)
            mtag_lepton_iter.append(1)
            mask_lepton_iter.append(True)
            phi_lepton_iter.append(event.mu_phi[i])
            pt_lepton_iter.append(event.mu_pt_NOSYS[i])
            e_lepton_iter.append(event.mu_e[i])
            q_lepton_iter.append(event.mu_charge[i])
            
        # No do our jets   
        for j in range(maxjets):
            if j >= len(event.jet_eta):
                mask_jet_iter.append(False)
                etai.append(0)
                btagi.append(0)
                phii.append(0)
                pti.append(0)
                ei.append(0)
            else:
                mask_jet_iter.append(True)
                eta_jet_iter.append(event.jet_eta[j])
                btag_jet_iter.append(int(event.jet_GN2v01_FixedCutBEff_85_select[j]==chr(1)))
                phi_jet_iter.append(event.jet_phi[j])
                pt_jet_iter.append(event.jet_pt_NOSYS[j])
                e_jet_iter.append(event.jet_e[j])


        mask_jet.append(mask_jet_iter)
        eta_jet.append(eta_jet_iter)
        btag_jet.append(btag_jet_iter)
        phi_jet.append(phi_jet_iter)
        pt_jet.append(pt_jet_iter)
        e_jet.append(e_jet_iter)

        mask_lepton.append(mask_lepton_iter)
        eta_lepton.append(eta_lepton_iter)
        etag_lepton.append(etag_lepton_iter)
        mtag_lepton.append(mtag_lepton_iter)
        q_lepton.append(q_lepton_iter)
        phi_lepton.append(phi_lepton_iter)
        pt_lepton.append(pt_lepton_iter)
        e_lepton.append(e_lepton_iter)

        first_lt = 0
        first_ht = 0
        for t in range(4):
            if event.parton_top_isHadronic[t] == 1:
                if first_ht == 0:
                    b1.append(event.b_recoj_index[t] if event.b_recoj_index[t] != -1 else -1)
                    q11.append(event.wd1_recoj_index[t] if event.wd1_recoj_index[t] != -1 else -1)
                    q12.append(event.wd2_recoj_index[t] if event.wd2_recoj_index[t] != -1 else -1)

                    #check for jet matching errors
                    if event.b_recoj_index[t] == event.wd1_recoj_index[t] and event.wd1_recoj_index[t] != -1:
                        b1[-1] = -1
                        q11[-1] = -1
                    if event.wd1_recoj_index[t] == event.wd2_recoj_index[t] and event.wd1_recoj_index[t] != -1:
                        q11[-1] = -1
                        q12[-1] = -1
                    if event.b_recoj_index[t] == event.wd2_recoj_index[t] and event.wd2_recoj_index[t] != -1:
                        b1[-1] = -1
                        q12[-1] = -1
                    first_ht = 1
                else:
                    b2.append(event.b_recoj_index[t] if event.b_recoj_index[t] != -1 else -1)
                    q21.append(event.wd1_recoj_index[t] if event.wd1_recoj_index[t] != -1 else -1)
                    q22.append(event.wd2_recoj_index[t] if event.wd2_recoj_index[t] != -1 else -1)
                    
                    # check for matching errors
                    if event.b_recoj_index[t] == event.wd1_recoj_index[t] and event.wd1_recoj_index[t] != -1:
                        b2[-1] = -1
                        q21[-1] = -1
                    if event.wd1_recoj_index[t] == event.wd2_recoj_index[t] and event.wd1_recoj_index[t] != -1:
                        q21[-1] = -1
                        q22[-1] = -1
                    if event.b_recoj_index[t] == event.wd2_recoj_index[t] and event.wd2_recoj_index[t] != -1:
                        b2[-1] = -1
                        q22[-1] = -1
                        
            else:
                if first_lt == 0:
                    b3.append(event.b_recoj_index[t] if event.b_recoj_index[t] != -1 else -1)
        
                    if event.wd1_el_index[t] != -1:
                        l3.append(event.wd1_el_index[t])
                    elif event.wd2_el_index[t] != -1:
                        l3.append(event.wd2_el_index[t])
                    elif event.wd1_mu_index[t] != -1:
                        l3.append(event.wd1_mu_index[t] + event.nElectrons)
                    elif event.wd2_mu_index[t] != -1:
                        l3.append(event.wd2_mu_index[t] + event.nElectrons)
                    else:
                        l3.append(-1)
                    first_lt = 1
                else:
                    b4.append(event.b_recoj_index[t] if event.b_recoj_index[t] != -1 else -1)
                    if event.wd1_el_index[t] != -1:
                        l4.append(event.wd1_el_index[t])
                    elif event.wd2_el_index[t] != -1:
                        l4.append(event.wd2_el_index[t])
                    elif event.wd1_mu_index[t] != -1:
                        l4.append(event.wd1_mu_index[t] + event.nElectrons)
                    elif event.wd2_mu_index[t] != -1:
                        l4.append(event.wd2_mu_index[t] + event.nElectrons)
                    else:
                        l4.append(-1)

                    


jet = inputs.create_group('Jet')
mask_data_jet = jet.create_dataset('MASK', data=mask_jet)
eta_data_jet = jet.create_dataset('eta', data=eta_jet)
btag_data_jet = jet.create_dataset('btag', data=btag_jet)
phi_data_jet = jet.create_dataset('phi', data=phi_jet)
pt_data_jet = jet.create_dataset('pt', data=pt_jet)
e_data_jet = jet.create_dataset('e', data=e_jet)

lepton = inputs.create_group('Lepton')
mask_data_lepton = lepton.create_dataset('MASK', data=mask_lepton)
eta_data_lepton = lepton.create_dataset('eta', data=eta_lepton)
btag_data_lepton = lepton.create_dataset('btag', data=btag_lepton)
phi_data_lepton = lepton.create_dataset('phi', data=phi_lepton)
pt_data_lepton = lepton.create_dataset('pt', data=pt_lepton)
e_data_lepton = lepton.create_dataset('e', data=e_lepton)


b3_data = t3.create_dataset('b', data=b3)
l3_data = t3.create_dataset('l', data=l3)
b4_data = t4.create_dataset('b', data=b4)
l4_data = t4.create_dataset('l', data=l4)

b1_data = t1.create_dataset('b', data=b1)
q11_data = t1.create_dataset('q1', data=q11)
q12_data = t1.create_dataset('q2', data=q12)
b2_data = t2.create_dataset('b', data=b2)
q21_data = t2.create_dataset('q1', data=q21)
q22_data = t2.create_dataset('q2', data=q22)


hf.close()

#user.dreiter.38180875._000001.output.root



#





