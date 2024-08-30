import ROOT
import numpy as np
import sys

canvas = ROOT.TCanvas("canvas", "output", 1800, 1200)
file = ROOT.TFile.Open(str(sys.argv[1]))
tree = file.Get("reco")

mt_hist = ROOT.TH1D("mt_hist","SPANet Hadronic Top Quark Mass (tttt 2LSS All)",100,0,500000)
mW_hist = ROOT.TH1D("mW_hist","SPANet Hadronic W Boson Mass (tttt 2LSS All)",100,0,400000)
mbl_hist = ROOT.TH1D("mW_hist","SPANet B-Jet + Lepton Mass (tttt 2LSS All)",100,0,400000)

ap_hist = ROOT.TH1D("ap_hist","SPANet Assignment Probability (tttt 2LSS All)",50,10,10)
dp_hist = ROOT.TH1D("dp_hist","SPANet Detection Probability (tttt 2LSS All)",50,10,10)
mp_hist = ROOT.TH1D("mp_hist","SPANet Marginal Probability (tttt 2LSS All)",50,10,10)


m_prev = 0
for event in tree:
    if event.top_m_SPANET[0] == m_prev and m_prev != -1 and m_prev != 0:
        break
    else:
        m_prev = event.top_m_SPANET[0]
    
    for t in range(2):
        m = event.top_m_SPANET[t]
        if m != -1:
            mt_hist.Fill(m)
    for t in range(2):
        m = event.W_m_SPANET[t]
        if m != -1:
            mW_hist.Fill(m)
    for t in range(2,4):
        m = event.mbl_SPANET[t]
        if m != -1:
            mbl_hist.Fill(m)
    
    for p in event.top_assign_prob_SPANET:
        ap_hist.Fill(p)
    for p in event.top_detect_prob_SPANET:
        dp_hist.Fill(p)
    for p in event.top_margin_prob_SPANET:
        mp_hist.Fill(p)

mt_hist.GetXaxis().SetTitle("m_{t} [MeV]")
mt_hist.GetYaxis().SetTitle("Candidates")
mt_hist.Draw()
canvas.Draw()
canvas.SaveAs("tttt_2lss_multi_hists/mt_hist.png")
canvas.Clear()

mW_hist.GetXaxis().SetTitle("m_{W} [MeV]")
mW_hist.GetYaxis().SetTitle("Candidates")
mW_hist.Draw()
canvas.Draw()
canvas.SaveAs("tttt_2lss_multi_hists/mW_hist.png")
canvas.Clear()

mbl_hist.GetXaxis().SetTitle("m_{bl} [MeV]")
mbl_hist.GetYaxis().SetTitle("Candidates")
mbl_hist.Draw()
canvas.Draw()
canvas.SaveAs("tttt_2lss_multi_hists/mbl_hist.png")
canvas.Clear()

ap_hist.GetYaxis().SetTitle("Candidates")
ap_hist.Draw()
canvas.Draw()
canvas.SaveAs("tttt_2lss_multi_hists/ap_hist.png")
canvas.Clear()

dp_hist.GetYaxis().SetTitle("Candidates")
dp_hist.Draw()
canvas.Draw()
canvas.SaveAs("tttt_2lss_multi_hists/dp_hist.png")
canvas.Clear()

mp_hist.GetYaxis().SetTitle("Candidates")
mp_hist.Draw()
canvas.Draw()
canvas.SaveAs("tttt_2lss_multi_hists/mp_hist.png")
canvas.Clear()
