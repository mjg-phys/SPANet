import ROOT
from ROOT import Math
import h5py
import sys
import numpy as np
from array import array

canvas = ROOT.TCanvas("canvas", "output", 1800, 1200)

mt_hist = ROOT.TH1D("mt_hist","SPANet Hadronic Top Quark Mass (tttt 2LSS All)",100,0,500000)
mW_hist = ROOT.TH1D("mW_hist","SPANet Hadronic W Boson Mass (tttt 2LSS All)",100,0,400000)

ap_hist = ROOT.TH1D("ap_hist","SPANet Assignment Probability (tttt 2LSS All)",50,10,10)
dp_hist = ROOT.TH1D("dp_hist","SPANet Detection Probability (tttt 2LSS All)",50,10,10)
mp_hist = ROOT.TH1D("mp_hist","SPANet Marginal Probability (tttt 2LSS All)",50,10,10)

hf = h5py.File("/eos/user/d/dreiter/SPANet/predict/tttt_2lss_hadronic_testing_output2.h5", 'r')

output_folder = "/eos/user/d/dreiter/SPANet/predict/tttt_2lss_hadronic_hists/"

def t_m(pt, eta, phi, e, b_index, j1_index, j2_index):
    # print(b_index-2)
    # print(j1_index-2)
    # print(j2_index-2)
    # print("njets " + str(len(event.jet_pt_NOSYS)))
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

b1 = np.array(hf.get('TARGETS/t1/b'))
q11 = np.array(hf.get('TARGETS/t1/q1'))
q12 = np.array(hf.get('TARGETS/t1/q2'))

b2 = np.array(hf.get('TARGETS/t2/b'))
q21 = np.array(hf.get('TARGETS/t2/q1'))
q22 = np.array(hf.get('TARGETS/t2/q2'))

ap1 = np.array(hf.get('TARGETS/t1/assignment_probability'))
ap2 = np.array(hf.get('TARGETS/t2/assignment_probability'))

dp1 = np.array(hf.get('TARGETS/t1/detection_probability'))
dp2 = np.array(hf.get('TARGETS/t2/detection_probability'))

mp1 = np.array(hf.get('TARGETS/t1/marginal_probability'))
mp2 = np.array(hf.get('TARGETS/t2/marginal_probability'))

pt = np.array(hf.get('INPUTS/Source/pt'))
eta = np.array(hf.get('INPUTS/Source/eta'))
phi = np.array(hf.get('INPUTS/Source/phi'))
e = np.array(hf.get('INPUTS/Source/e'))


for i in range(len(b1)): 


    tm1 = t_m(pt[i], eta[i], phi[i], e[i], b1[i], q11[i], q12[i])
    tm2 = t_m(pt[i], eta[i], phi[i], e[i], b2[i], q21[i], q22[i])
    
    wm1 = W_m(pt[i], eta[i], phi[i], e[i], q11[i], q12[i])
    wm2 = W_m(pt[i], eta[i], phi[i], e[i], q21[i], q22[i])

    mt_hist.Fill(tm1)
    mt_hist.Fill(tm2)

    mW_hist.Fill(wm1)
    mW_hist.Fill(wm2)

    ap_hist.Fill(ap1[i])
    ap_hist.Fill(ap2[i])

    dp_hist.Fill(dp1[i])
    dp_hist.Fill(dp2[i])

    mp_hist.Fill(mp1[i])
    mp_hist.Fill(mp2[i])

mt_hist.GetXaxis().SetTitle("m_{t} [MeV]")
mt_hist.GetYaxis().SetTitle("Candidates")
mt_hist.Draw()
canvas.Draw()
canvas.SaveAs(output_folder + "mt_hist.png")
canvas.Clear()

mW_hist.GetXaxis().SetTitle("m_{W} [MeV]")
mW_hist.GetYaxis().SetTitle("Candidates")
mW_hist.Draw()
canvas.Draw()
canvas.SaveAs(output_folder + "mW_hist.png")
canvas.Clear()

ap_hist.GetYaxis().SetTitle("Candidates")
ap_hist.Draw()
canvas.Draw()
canvas.SaveAs(output_folder + "ap_hist.png")
canvas.Clear()

dp_hist.GetYaxis().SetTitle("Candidates")
dp_hist.Draw()
canvas.Draw()
canvas.SaveAs(output_folder + "dp_hist.png")
canvas.Clear()

mp_hist.GetYaxis().SetTitle("Candidates")
mp_hist.Draw()
canvas.Draw()
canvas.SaveAs(output_folder + "mp_hist.png")
canvas.Clear()