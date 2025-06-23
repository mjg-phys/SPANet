import ROOT
from ROOT import Math
import h5py
import sys
import numpy as np
from array import array
import glob

def t_m(pt, eta, phi, e, b_index, j1_index, j2_index):
    b = Math.PtEtaPhiEVector()
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
    
    # if (bl.M() == 0):
    #     return np.nan
    
    return bl.M()


hf = h5py.File(str(sys.argv[1]), 'r+')

b1 = np.array(hf.get('TARGETS/t1/b'))
q11 = np.array(hf.get('TARGETS/t1/q1'))
q12 = np.array(hf.get('TARGETS/t1/q2'))

b2 = np.array(hf.get('TARGETS/t2/b'))
q21 = np.array(hf.get('TARGETS/t2/q1'))
q22 = np.array(hf.get('TARGETS/t2/q2'))

b3 = np.array(hf.get('TARGETS/t3/b'))
l3 = np.array(hf.get('TARGETS/t3/l'))

b4 = np.array(hf.get('TARGETS/t4/b'))
l4 = np.array(hf.get('TARGETS/t4/l'))

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
lepton_pt = np.pad(lepton_pt, ((0, 0), (16, 0)), mode='constant', constant_values=0)
lepton_eta = np.pad(lepton_eta, ((0, 0), (16, 0)), mode='constant', constant_values=0)
lepton_phi = np.pad(lepton_phi, ((0, 0), (16, 0)), mode='constant', constant_values=0)
lepton_e = np.pad(lepton_e, ((0, 0), (16, 0)), mode='constant', constant_values=0)


b1_og = np.array(hf.get('TARGETS_ORIGINAL/t1/b'))

# regressions_spanet = hf.create_group('REGRESSIONS_SPANET')

# t1 = regressions_spanet.create_group('t1')
# t2 = regressions_spanet.create_group('t2')
# t3 = regressions_spanet.create_group('t3')
# t4 = regressions_spanet.create_group('t4')


print(len(lepton_pt))

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


top_m_SPANET = []
top_isHadronic_SPANET = []
W_m_SPANET = []
mbl_SPANET = []

for i in range(len(b1)): 
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
        if (q12[i]>15):
            q12[i] =15
        if (q11[i]>15):
            q11[i] =15  
                    
            # continue
        b_indices = [b1[i], b2[i], b3[i], b4[i]]
        # print("b_indices: ", b_indices)

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
        
    top_m_SPANET.append(top_m)
    W_m_SPANET.append(W_mass)
    mbl_SPANET.append(m_bl)

    # for m in top_m:
    #     top_m_SPANET.push_back(m)
        
    # for m in W_mass:
    #     W_m_SPANET.push_back(m)

    # for t in top_had:
    #     top_isHadronic_SPANET.push_back(t)

    # for j in m_bl:
    #     mbl_SPANET.push_back(j)

print(hf)
regression_group = hf['REGRESSIONS']
t1_reg = regression_group["t1"]
t2_reg = regression_group["t2"]
t3_reg = regression_group["t3"]
t4_reg = regression_group["t4"]

t1_mass = [arr[0] for arr in top_m_SPANET]
t2_mass = [arr[1] for arr in top_m_SPANET]

w1_mass = [arr[0] for arr in W_m_SPANET]
w2_mass = [arr[1] for arr in W_m_SPANET]

mbl3 = [arr[2] for arr in mbl_SPANET]
mbl4 = [arr[3] for arr in mbl_SPANET]

t1_reg.create_dataset('tMass_SPANET', data=t1_mass)
t2_reg.create_dataset('tMass_SPANET', data=t2_mass)

t1_reg.create_dataset('wMass_SPANET', data=w1_mass)
t2_reg.create_dataset('wMass_SPANET', data=w2_mass)

t3_reg.create_dataset('mbl_SPANET', data=mbl3)
print(mbl3)
t4_reg.create_dataset('mbl_SPANET', data=mbl4)

print(mbl4)

