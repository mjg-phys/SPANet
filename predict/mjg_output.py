import ROOT
from ROOT import Math
import h5py
import sys
import numpy as np
from array import array

hf = h5py.File(str(sys.argv[2]), 'r')

# chain = ROOT.TChain("reco")
# chain.Add(str(sys.argv[1]))
# tree = chain.CopyTree("(pass_SSee_passECIDS_NOSYS||pass_SSem_passECIDS_NOSYS||pass_SSmm_NOSYS)&&(!(pass_eee_ZVeto_NOSYS||pass_eem_ZVeto_NOSYS||pass_emm_ZVeto_NOSYS||pass_mmm_ZVeto_NOSYS))&&(!pass_llll_NOSYS)")

# output_file = ROOT.TFile.Open(str(sys.argv[3]), "RECREATE")



chain = ROOT.TChain("reco")
chain.Add(str(sys.argv[1]))
tree = chain.CopyTree("(pass_SSee_passECIDS_NOSYS||pass_SSem_passECIDS_NOSYS||pass_SSmm_NOSYS)&&(!(pass_eee_ZVeto_NOSYS||pass_eem_ZVeto_NOSYS||pass_emm_ZVeto_NOSYS||pass_mmm_ZVeto_NOSYS))&&(!pass_llll_NOSYS)")
new_tree = tree.CloneTree(0)  # Create an empty clone with the same structure

output_file = ROOT.TFile.Open(str(sys.argv[3]), "RECREATE")
# tree = output_file.Get("reco")
# num_entries = tree.GetEntries()
# print(f"The tree contains {num_entries} events.")
# exit()


top_m_SPANET = ROOT.std.vector[float]()
top_isHadronic_SPANET = ROOT.std.vector[int]()
W_m_SPANET = ROOT.std.vector[float]()

b_index_SPANET = ROOT.std.vector[int]()
q1_index_SPANET = ROOT.std.vector[int]()
q2_index_SPANET = ROOT.std.vector[int]()

mbl_SPANET = ROOT.std.vector[float]()
el_index_SPANET = ROOT.std.vector[int]()
mu_index_SPANET = ROOT.std.vector[int]()

top_assign_prob_SPANET = ROOT.std.vector[float]()
top_detect_prob_SPANET = ROOT.std.vector[float]()
top_margin_prob_SPANET = ROOT.std.vector[float]()

new_tree.Branch("top_m_SPANET", top_m_SPANET)
new_tree.Branch("top_isHadronic_SPANET", top_isHadronic_SPANET)
new_tree.Branch("W_m_SPANET", W_m_SPANET)

new_tree.Branch("b_index_SPANET",  b_index_SPANET)
new_tree.Branch("q1_index_SPANET",  q1_index_SPANET)
new_tree.Branch("q2_index_SPANET", q2_index_SPANET)

new_tree.Branch("mbl_SPANET",  mbl_SPANET)
new_tree.Branch("el_index_SPANET", el_index_SPANET)
new_tree.Branch("mu_index_SPANET",  mu_index_SPANET)

new_tree.Branch("top_assign_prob_SPANET", top_assign_prob_SPANET)
new_tree.Branch("top_detect_prob_SPANET",  top_detect_prob_SPANET)
new_tree.Branch("top_margin_prob_SPANET", top_margin_prob_SPANET)

def t_m(pt, eta, phi, e, b_index, j1_index, j2_index):
    b = Math.PtEtaPhiEVector()
    print(b_index)
    b.SetCoordinates(pt[b_index], eta[b_index], phi[b_index], e[b_index])

    j1 = Math.PtEtaPhiEVector()
    j1.SetCoordinates(pt[j1_index], eta[j1_index], phi[j1_index], e[j1_index])

    j2 = Math.PtEtaPhiEVector()
    j2.SetCoordinates(pt[j2_index], eta[j2_index], phi[j2_index], e[j2_index])
            
    t = b + j1 + j2
    return t.M()

def W_m(pt, eta, phi, e, j1_index, j2_index):
    j1 = Math.PtEtaPhiEVector()
    j1.SetCoordinates(pt[j1_index], eta[j1_index], phi[j1_index], e[j1_index])

    j2 = Math.PtEtaPhiEVector()
    j2.SetCoordinates(pt[j2_index], eta[j2_index], phi[j2_index], e[j2_index])

    W = j1 + j2
    return W.M()

def mbl(pt, eta, phi, e, lepton_pt,lepton_eta,lepton_phi,lepton_e, b_index, l_index):
    b = Math.PtEtaPhiEVector()
    b.SetCoordinates(pt[b_index], eta[b_index], phi[b_index], e[b_index])

   
    l = Math.PtEtaPhiEVector()
    l.SetCoordinates(lepton_pt[l_index], lepton_eta[l_index], lepton_phi[l_index],lepton_e[l_index])
    
    bl = b + l
    return bl.M()

i = 0

b1 = np.array(hf.get('TARGETS/t1/b'))
q11 = np.array(hf.get('TARGETS/t1/q1'))
q12 = np.array(hf.get('TARGETS/t1/q2'))

b2 = np.array(hf.get('TARGETS/t2/b'))
q21 = np.array(hf.get('TARGETS/t2/q1'))
q22 = np.array(hf.get('TARGETS/t2/q2'))

b3 = np.array(hf.get('TARGETS/t3/b'))
print(b3)
l3 = np.array(hf.get('TARGETS/t3/l'))
print(l3)
print(min(l3))

b4 = np.array(hf.get('TARGETS/t4/b'))
l4 = np.array(hf.get('TARGETS/t4/l'))
print(min(l4))
print(l3)
print(l4)
ap = np.array([np.array(hf.get('TARGETS/t1/assignment_probability')), np.array(hf.get('TARGETS/t2/assignment_probability')), np.array(hf.get('TARGETS/t3/assignment_probability')), np.array(hf.get('TARGETS/t4/assignment_probability'))])
dp = np.array([np.array(hf.get('TARGETS/t1/detection_probability')), np.array(hf.get('TARGETS/t2/detection_probability')), np.array(hf.get('TARGETS/t3/detection_probability')), np.array(hf.get('TARGETS/t4/detection_probability'))])
mp = np.array([np.array(hf.get('TARGETS/t1/marginal_probability')), np.array(hf.get('TARGETS/t2/marginal_probability')), np.array(hf.get('TARGETS/t3/marginal_probability')), np.array(hf.get('TARGETS/t4/marginal_probability'))])

pt = np.array(hf.get('INPUTS/Jet/pt'))
eta = np.array(hf.get('INPUTS/Jet/eta'))
phi = np.array(hf.get('INPUTS/Jet/phi'))
e = np.array(hf.get('INPUTS/Jet/e'))

lepton_pt = np.array(hf.get('INPUTS/Lepton/pt'))
lepton_eta = np.array(hf.get('INPUTS/Lepton/eta'))
lepton_phi = np.array(hf.get('INPUTS/Lepton/phi'))
lepton_e = np.array(hf.get('INPUTS/Lepton/e'))
# l3 = l3-16
# l4 = l4-16

print(len(lepton_pt))

# For some reason some events have gone rogue and set 
for j in range(len(l3)):
    if (l3[j] ==16 or l4[j]==17):
        l3[j] = 1
        l4[j] = 0
    else:
        l3[j] = 0
        l4[j] =1

# print(max(l3))
# exit()

for event in tree: 
    b_indices = [-1,-1,-1,-1]
    q1_indices = [-1,-1,-1,-1]
    q2_indices = [-1,-1,-1,-1]
    top_m = [-1,-1,-1,-1]
    W_mass = [-1,-1,-1,-1]
    top_had = [-1,-1,-1,-1]
    mu_index = [-1,-1,-1,-1]
    el_index = [-1,-1,-1,-1]
    top_assign_prob = [-1,-1,-1,-1]
    top_detect_prob = [-1,-1,-1,-1]
    top_margin_prob = [-1,-1,-1,-1]
    if (i >= len(b1)):    
        break
    else:
        if (b1[i]>15):
            b1[i] =15
        if (b2[i]>15):
            b2[i] =15
        if (b3[i]>15):
            b3[i] =15
        if (b4[i]>15):
            b4[i] =15
        if (q22[i]>15):
            q22[i] =15
        if (q21[i]>15):
            q21[i] =15           
            # continue
        print(pt)
        print(i)
        print( b1[i])
        print(b2[i])
        
        tm1 = t_m(pt[i], eta[i], phi[i], e[i], b1[i], q11[i], q12[i])
        tm2 = t_m(pt[i], eta[i], phi[i], e[i], b2[i], q21[i], q22[i])
        # top_m = [tm1 if b1[i] >= 2 and q11[i] >= 2 and q12[i] >= 2 and tm1 > 0 else -1, tm2 if b2[i] >= 2 and q21[i] >= 2 and q22[i] >= 2 and tm2 > 0 else -1]
        top_m = [tm1,tm2,-1,-1]

        wm1 = W_m(pt[i], eta[i], phi[i], e[i], q11[i], q12[i])
        wm2 = W_m(pt[i], eta[i], phi[i], e[i], q21[i], q22[i])
        W_mass = [wm1, wm2 , -1, -1]
        top_had = [1,1,0,0]

        etag = np.array(hf.get('INPUTS/Source/etag'))
        mtag = np.array(hf.get('INPUTS/Source/mtag'))

        mbl3 = mbl(pt[i], eta[i], phi[i], e[i], lepton_pt[i], lepton_eta[i], lepton_phi[i], lepton_e[i], b3[i], l3[i])
        mbl4 = mbl(pt[i], eta[i], phi[i], e[i], lepton_pt[i], lepton_eta[i], lepton_phi[i], lepton_e[i], b4[i], l4[i])
        m_bl = [-1, -1, mbl3, mbl4]

        # el_index = [-1, -1, int(l3[i]) if etag[i][l3[i]] == 1 and l3[i] < 2 else -1, int(l4[i]) if etag[i][l4[i]] == 1 and l4[i] < 2 else -1]
        
        # nElectrons = etag[i][0] + etag[i][1]
        # nMuons = mtag[i][0] + mtag[i][1]

        # mu_index = [-1, -1, int(l3[i]-nElectrons) if mtag[i][l3[i]] == 1 and l3[i] < 2 else -1, int(l4[i]-nElectrons) if mtag[i][l4[i]] == 1 and l4[i] < 2 else -1]

        top_assign_prob = [ap[t][i] for t in range(4)]

        top_detect_prob = [dp[t][i] for t in range(4)]

        top_margin_prob = [mp[t][i] for t in range(4)]


    top_m_SPANET.clear()
    top_m_SPANET.reserve(4)
    for m in top_m:
        top_m_SPANET.push_back(m)
        
    W_m_SPANET.clear()
    W_m_SPANET.reserve(4)
    for m in W_mass:
        W_m_SPANET.push_back(m)

    top_isHadronic_SPANET.clear()
    top_isHadronic_SPANET.reserve(4)
    for t in top_had:
        top_isHadronic_SPANET.push_back(t)

    b_index_SPANET.clear()
    b_index_SPANET.reserve(4)
    for j in b_indices:
        b_index_SPANET.push_back(j)

    q1_index_SPANET.clear()
    q1_index_SPANET.reserve(4)
    for j in q1_indices:
        q1_index_SPANET.push_back(j)

    q2_index_SPANET.clear()
    q2_index_SPANET.reserve(4)
    for j in q2_indices:
        q2_index_SPANET.push_back(j)

    mbl_SPANET.clear()
    mbl_SPANET.reserve(4)
    for j in m_bl:
        mbl_SPANET.push_back(j)

    # el_index_SPANET.clear()
    # el_index_SPANET.reserve(4)
    # for j in el_index:
    #     el_index_SPANET.push_back(j)    
    
    # mu_index_SPANET.clear()
    # mu_index_SPANET.reserve(4)
    # for j in mu_index:
    #     mu_index_SPANET.push_back(j)  

    top_assign_prob_SPANET.clear()
    top_assign_prob_SPANET.reserve(4)
    for j in top_assign_prob:
        top_assign_prob_SPANET.push_back(j)

    top_detect_prob_SPANET.clear()
    top_detect_prob_SPANET.reserve(4)
    for j in top_detect_prob:
        top_detect_prob_SPANET.push_back(j)

    top_margin_prob_SPANET.clear()
    top_margin_prob_SPANET.reserve(4)
    for j in top_margin_prob:
        top_margin_prob_SPANET.push_back(j)

    # print(top_m)
    new_tree.Fill()      
    i+=1

new_tree.Write()
output_file.Close()