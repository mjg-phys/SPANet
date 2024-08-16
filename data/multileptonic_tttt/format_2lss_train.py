import h5py
import sys
import ROOT
import numpy as np

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

for n in sys.argv[2:]:
    chain.Add(str(n))

tree = chain.CopyTree("pass_SSee_passECIDS_NOSYS||pass_SSem_passECIDS_NOSYS||pass_SSmm_NOSYS")

name = sys.argv[1]

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
b3, l3, b4, l4 = [], [], [], []
b1, q11, q12, b2, q21, q22 = [], [], [], [], [], []

for event in tree:

    if len(event.parton_top_isHadronic < 4) or sum(event.parton_top_isHadronic) != 2:
        continue
    
    if sum(event.parton_top_isHadronic) == 2 and len(event.parton_top_isHadronic) == 4:
        maski, etai, etagi, mtagi, btagi, phii, pti, ei, qi = [], [], [],[], [], [], [], [], []
        for i in range(event.nElectrons):
            etai.append(event.el_eta[i])
            etagi.append(1)
            mtagi.append(0)
            btagi.append(0)
            maski.append(True)
            phii.append(event.el_phi[i])
            pti.append(event.el_pt_NOSYS[i])
            ei.append(event.el_e[i])
            qi.append(event.el_charge[i])
        for i in range(event.nMuons):
            etai.append(event.mu_eta[i])
            etagi.append(0)
            mtagi.append(1)
            btagi.append(0)
            maski.append(True)
            phii.append(event.mu_phi[i])
            pti.append(event.mu_pt_NOSYS[i])
            ei.append(event.mu_e[i])
            qi.append(event.mu_charge[i])
        for j in range(maxjets):
            qi.append(0)
            mtagi.append(0)
            etagi.append(0)
            if j >= len(event.jet_eta):
                maski.append(False)
                etai.append(0)
                btagi.append(0)
                phii.append(0)
                pti.append(0)
                ei.append(0)
            else:
                maski.append(True)
                etai.append(event.jet_eta[j])
                btagi.append(int(event.jet_GN2v01_FixedCutBEff_85_select[j]==chr(1)))
                phii.append(event.jet_phi[j])
                pti.append(event.jet_pt_NOSYS[j])
                ei.append(event.jet_e[j])

        eta.append(etai)
        etag.append(etagi)
        mtag.append(mtagi)
        mask.append(maski)
        btag.append(btagi)
        phi.append(phii)
        pt.append(pti)
        e.append(ei)
        q.append(qi)

        first_lt = 0
        first_ht = 0
        for t in range(4):
            if event.parton_top_isHadronic[t] == 1:
                if first_ht == 0:
                    b1.append(event.b_recoj_index[t]+2 if event.b_recoj_index[t] != -1 else -1)
                    q11.append(event.wd1_recoj_index[t]+2 if event.wd1_recoj_index[t] != -1 else -1)
                    q12.append(event.wd2_recoj_index[t]+2 if event.wd2_recoj_index[t] != -1 else -1)

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
                    b2.append(event.b_recoj_index[t]+2 if event.b_recoj_index[t] != -1 else -1)
                    q21.append(event.wd1_recoj_index[t]+2 if event.wd1_recoj_index[t] != -1 else -1)
                    q22.append(event.wd2_recoj_index[t]+2 if event.wd2_recoj_index[t] != -1 else -1)
                    
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
                    b3.append(event.b_recoj_index[t]+2 if event.b_recoj_index[t] != -1 else -1)
        
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
                    b4.append(event.b_recoj_index[t]+2 if event.b_recoj_index[t] != -1 else -1)
                    
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

                    


source = inputs.create_group('Source')
mask_data = source.create_dataset('MASK', data=mask)
eta_data = source.create_dataset('eta', data=eta)
btag_data = source.create_dataset('btag', data=btag)
phi_data = source.create_dataset('phi', data=phi)
pt_data = source.create_dataset('pt', data=pt)
e_data = source.create_dataset('e', data=e)
q_data = source.create_dataset('q', data=q)
etag_data = source.create_dataset('etag', data=etag)
mtag_data = source.create_dataset('mtag', data=mtag)


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





