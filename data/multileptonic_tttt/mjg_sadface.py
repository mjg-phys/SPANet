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

# for n in sys.argv[2:]:
#     chain.Add(str(n))
# Find all .root files in the folder

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

# Fill data lists
mask, eta, etag, btag, mtag, phi, pt, e, q = [], [], [],[], [], [], [], [], []
jet_mask, jet_eta, jet_etag, jet_btag, jet_mtag, jet_phi, jet_pt, jet_e, jet_q = [], [], [],[], [], [], [], [], []
lepton_mask, lepton_eta, lepton_etag, lepton_btag, lepton_mtag, lepton_phi, lepton_pt, lepton_e, lepton_q = [], [], [],[], [], [], [], [], []


b3, l3, b4, l4 = [], [], [], []
b1, q11, q12, b2, q21, q22 = [], [], [], [], [], []

for event in tree:
    if len(event.parton_top_isHadronic) != 4:
        continue
    if sum(event.parton_top_isHadronic) != 2:
        continue
    
    lepton_maski, lepton_etai, lepton_etagi, lepton_mtagi, lepton_btagi, lepton_phii, lepton_pti, lepton_ei, lepton_qi = [], [], [],[], [], [], [], [], []
    
    jet_maski, jet_etai, jet_etagi, jet_mtagi, jet_btagi, jet_phii, jet_pti, jet_ei, jet_qi = [], [], [],[], [], [], [], [], []

    for i in range(event.nElectrons):
        lepton_etai.append(event.el_eta[i])
        lepton_etagi.append(1)
        lepton_mtagi.append(0)
        lepton_btagi.append(0)
        lepton_maski.append(True)
        lepton_phii.append(event.el_phi[i])
        lepton_pti.append(event.el_pt_NOSYS[i])
        lepton_ei.append(event.el_e[i])
        lepton_qi.append(event.el_charge[i])
    for i in range(event.nMuons):
        lepton_etai.append(event.mu_eta[i])
        lepton_etagi.append(0)
        lepton_mtagi.append(1)
        lepton_btagi.append(0)
        lepton_maski.append(True)
        lepton_phii.append(event.mu_phi[i])
        lepton_pti.append(event.mu_pt_NOSYS[i])
        lepton_ei.append(event.mu_e[i])
        lepton_qi.append(event.mu_charge[i])
    for j in range(maxjets):
        if j >= len(event.jet_eta):
            jet_maski.append(False)
            jet_etai.append(0)
            jet_btagi.append(0)
            jet_phii.append(0)
            jet_pti.append(0)
            jet_ei.append(0)
        else:
            jet_maski.append(True)
            jet_etai.append(event.jet_eta[j])
            jet_btagi.append(int(event.jet_GN2v01_FixedCutBEff_85_select[j]==chr(1)))
            jet_phii.append(event.jet_phi[j])
            jet_pti.append(event.jet_pt_NOSYS[j])
            jet_ei.append(event.jet_e[j])

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
    
    first_lt = 0
    first_ht = 0
    for t in range(4):
        # print(event.b_recoj_index)
        # print( event.parton_top_isHadronic[t])
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
    if len(b1) != len(b2):
        print(len(b1))
        print(len(b2))
        print(sum(event.parton_top_isHadronic))
        
        print("wtf" )
        exit()

                    
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

hf.close()





