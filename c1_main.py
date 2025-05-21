# %%
import os, sys, math, ROOT
from array import array
import numpy as np
from numpy import cos, sin, tan, cosh, sinh
from numpy.linalg import norm
from numpy.ma.core import arccos

from helper import *

#root_file = "data.root"
root_file = "GluGluHToZZTo4L_M125_new.root"

if root_file == "data.root":
    outputdirect = 'SIGNAL/FF'
elif root_file == "GluGluHToZZTo4L_M125_new.root":
    outputdirect = 'Signal/F4_TightID'
else:
    print("Invalid root_file --> no outputdirect")

# Input file
inputfile = ROOT.TFile(root_file, "READ")

# Output file
outputfile = ROOT.TFile("output_4m.root", "RECREATE")

# Accessing the "TTree" from root that contains the data (one entry per event)
tree = inputfile.Get("ntuplizer/tree")

nentries = 5000  # Run only on partial events to speed up the code
nentries = tree.GetEntries()
print("Number of entries: ", nentries)

hist_mass4m = ROOT.TH1F("hist_M_4muons", ";Invariant Mass of the four muons (GeV);Events", 125, 0, 500)
hist_N_muons = ROOT.TH1F("hist_N_of_muons", ";Number of muons;Events", 70, 0, 7)

hist_mass4m_z = ROOT.TH1F("hist_M_4muons", ";Invariant Mass of the four muons (GeV);Events", 50, 100, 150)
hist_mass4m_bkg = ROOT.TH1F("hist_background", ";Invariant Mass of the four muons (GeV);Events", 50, 100, 150)
hist_mass4m_pic = ROOT.TH1F("hist_higgs_peak", ";Invariant Mass of the four muons (GeV);Events", 50, 100, 150)

hist_mass2m_a = ROOT.TH1F("hist_m_pair1", ";Invariant Mass of muon pair a (GeV);Events", 40, 0, 20)
hist_mass2m_b = ROOT.TH1F("hist_m_pair2", ";Invariant Mass of muon pair b (GeV);Events", 40, 0, 20)
hist_mass2m_2m = ROOT.TH2F("hist_m_pairs", ";Invariant Mass of muon pair a (GeV);Invariant Mass of muon pair b (GeV);Events",
                         100, 0, 200, 100, 0, 200)
hist_mass4m_2m_a = ROOT.TH2F("hist_massmuonpair1&mass4muon", ";4 muon Invariant Mass;Invariant Mass of muon pair a;Events",
                         100, 0, 500, 100, 0, 200)
hist_mass4m_2m_b = ROOT.TH2F("hist_massmuonpair2&mass4muon", ";4 muon Invariant Mass;Invariant Mass of muon pair b;Events",
                         100, 0, 500, 100, 0, 200)
N_left = 0
N_right = 0
more_events = 0
less_events = 0
exact_events = 0


for i in range(0, nentries):
    tree.GetEntry(i)
    ''',
    TightID = 0
    for lep in range(0, len(tree._lPassTightID)):
        if tree._lPassTightID[lep] == 0:
            TightID += 1

    if TightID > 0:
        continue
    '''
    '''
    pdgId_absum = 0
    pdgId_sum = 0
    for lep in range(0, len(tree._lPt)):
        pdgId_sum += tree._lpdgId[lep]
        pdgId_absum += abs(tree._lpdgId[lep])

    # first selector skips the entry if not 2 muon-antimuon pairs
    if pdgId_sum != 0 or pdgId_absum != 52:
        continue
    '''

    muon = []
    antimuon = []
    for lep in range(0, len(tree._lPt)):
        #Check if the lepton indexed lep is a muon. If it is, add index in muon list
        if tree._lpdgId[lep] == 13:
            muon.append(lep)
        #Check if the lepton indexed lep is an antimuon. If it is, add index to antimuon list
        elif tree._lpdgId[lep] == -13:
            antimuon.append(lep)

    hist_N_muons.Fill(len(muon)+len(antimuon))

    if len(muon) < 2 or len(antimuon) < 2:
        less_events += 1

    if (len(muon) == 2) and (len(antimuon) == 2):
        exact_events += 1

    if (len(muon) != 2) or (len(antimuon) != 2):
        more_events += 1
        continue

    dist00 = np.sqrt((tree._lPhi[muon[0]]-tree._lPhi[antimuon[0]])**2
                        + (tree._lEta[muon[0]]-tree._lEta[antimuon[0]])**2)
    dist11 = np.sqrt((tree._lPhi[muon[1]]-tree._lPhi[antimuon[1]])**2
                        + (tree._lEta[muon[1]]-tree._lEta[antimuon[1]])**2)

    dist01 = np.sqrt((tree._lPhi[muon[0]]-tree._lPhi[antimuon[1]])**2
                        + (tree._lEta[muon[0]]-tree._lEta[antimuon[1]])**2)
    dist10 = np.sqrt((tree._lPhi[muon[1]]-tree._lPhi[antimuon[0]])**2
                        + (tree._lEta[muon[1]]-tree._lEta[antimuon[0]])**2)


    if (dist00 < dist01 and dist11 < dist10):
        pair1 = [muon[0], antimuon[0]]
        pair2 = [muon[1], antimuon[1]]
    elif dist00 + dist11 < dist01 + dist10:
        pair1 = [muon[0], antimuon[0]]
        pair2 = [muon[1], antimuon[1]]
    else:
        pair1 = [muon[0], antimuon[1]]
        pair2 = [muon[1], antimuon[0]]


    mass2m_1 = invariantmass2l(tree, pair1[0], pair1[1])
    mass2m_2 = invariantmass2l(tree, pair2[0], pair2[1])
    mass4m = invariantmass4l(tree, 0, 1, 2,3)

    # rearange pairs so that pair 'a' is always the most massive (~ 90GeV Z boson)
    mass2m_a = mass2m_1
    mass2m_b = mass2m_2

    if mass2m_1 > mass2m_2:
        pair_a = pair1
        pair_b = pair2
        mass2m_a = mass2m_1
        mass2m_b = mass2m_2
    else:
        pair_a = pair2
        pair_b = pair1
        mass2m_a = mass2m_2
        mass2m_b = mass2m_1

    res_ppZZ = 10

    # S4GeV
    if mass2m_a < 4 or mass2m_b < 4:
        continue

    # ppZZ
    if abs(mass2m_a-91) < res_ppZZ and abs(mass2m_b - 91) < res_ppZZ:
        continue

    # D20GeV
    #if mass2m_a+mass2m_b < 20 :
    #    continue

    # sorting the Pt values of the four leptons in ascending order
    pts = []
    pts_maxtomin = []

    for lep in range(0, len(tree._lPt)):
        pts.append(tree._lPt[lep])
    for lep in range(0, len(tree._lPt)):
        pts_maxtomin.append(max(pts))
        pts.remove(max(pts))

    # F3_Pt3_Pt4 : compare Pt3/4 of data vs simulation and cut accordingly
    if pts_maxtomin[2] < 15 or pts_maxtomin[3] < 10:
        continue
    
    SB_left = [104, 119]
    SB_right = [131, 146]
    res_Hw = 3
    if 100 < mass4m < 150:
        hist_mass4m_z.Fill(mass4m)
        if 125 - res_Hw < mass4m < 125 + res_Hw:
            hist_mass4m_pic.Fill(mass4m)
        if (SB_left[0] < mass4m < SB_left[1]) or (SB_right[0] < mass4m < SB_right[1]):
            hist_mass4m_bkg.Fill(mass4m)


    #hist_mass2m_a.Fill(mass2m_a)
    #hist_mass2m_b.Fill(mass2m_b)
    hist_mass4m.Fill(mass4m)
    #hist_mass2m_2m.Fill(mass2m_a, mass2m_b)
    #hist_mass4m_2m_a.Fill(mass4m, mass2m_a)
    #hist_mass4m_2m_b.Fill(mass4m, mass2m_b)



bkg_func = ROOT.TF1("bkg_func", "pol1", 105, 145)

hist_mass4m_bkg.Fit("bkg_func", "S")

fit = hist_mass4m_bkg.Fit("bkg_func", "S")
p0 = fit.Parameter(0)
p1 = fit.Parameter(1)
sig0 = fit.ParError(0)
sig1 = fit.ParError(1)
rho = fit.Correlation(0, 1)
print(rho)
'''
c = ROOT.TCanvas("hist_mass4m_pic", "hist_Higgs_peak", 600, 600)

hist_mass4m_pic.Draw()
bkg_func.Draw("SAME")
c.SaveAs(outputdirect + "/" + "hist_pic_and_bkg.pdf")
outputfile.cd()
c.Write()
c.Close()
'''
#savehisto(outputfile, hist_mass2m_a, "hist_mass2m_a_0-20GeV", outputdir=outputdirect)
#savehisto(outputfile, hist_mass2m_b, "hist_mass2m_b_x20", outputdir=outputdirect)
savehisto(outputfile, hist_mass4m, "hist_mass4m_FF", outputdir=outputdirect)
#savehisto(outputfile, hist_mass4m_z, "hist_mass4m_z", outputdir=outputdirect)
#savehisto(outputfile, hist_mass4m_bkg, "hist_bkg_04_18", outputdir=outputdirect)
#savehisto(outputfile, hist_mass4m_pic, "hist_pic_and_bkg", outputdir=outputdirect)
#savehisto(outputfile, hist_mass2m_2m, "hist2D_mass2m_2m", outputdir=outputdirect)
#savehisto(outputfile, hist_mass4m_2m_a, "hist2D_mass4m_2m_a", outputdir=outputdirect)
#savehisto(outputfile, hist_mass4m_2m_b, "hist2D_mass4m_2m_b", outputdir=outputdirect)
#savehisto(outputfile, hist_N_muons, 'hist_N_antimumu', outputdir=outputdirect)

#print(less_events, exact_events, more_events, less_events + exact_events + more_events)


