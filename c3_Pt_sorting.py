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
    outputdirect = '4MUON/Pt'
elif root_file == "GluGluHToZZTo4L_M125_new.root":
    outputdirect = 'SIGNAL/Pt_signal'
else:
    print("Invalid root_file => No outputdirect")

# Input file
inputfile = ROOT.TFile(root_file, "READ")

# Output file
outputfile = ROOT.TFile("output.root", "RECREATE")

# Accessing the "TTree" from root that contains the data (one entry per event)
tree = inputfile.Get("ntuplizer/tree")

nentries = 5000  # Run only on partial events to speed up the code
nentries = tree.GetEntries()
print("Number of entries: ", nentries)

hist_Pt1 = ROOT.TH1F("hist_Pt_1st", ";Highest Pt (GeV);Events", 50, 0, 100)
hist_Pt2 = ROOT.TH1F("hist_Pt_2nd", ";Second highest Pt (GeV);Events", 50, 0, 100)
hist_Pt3 = ROOT.TH1F("hist_Pt_3rd", ";Third highest Pt (GeV);Events", 50, 0, 100)
hist_Pt4 = ROOT.TH1F("hist_Pt_4th", ";Fourth highest Pt (GeV);Events", 50, 0, 100)

hist2D_Pt_am_1 = ROOT.TH2F("hist_Pt_muon_antimuon_1st",
                             ";Pt;Muon 0 - Antimuon 1;Events",
                         2, 0, 1, 50, 0, 50)
hist2D_Pt_am_2 = ROOT.TH2F("hist_Pt_muon_antimuon_2nd",
                             ";Pt;Muon 0 - Antimuon 1;Events",
                         2, 0, 1, 50, 0, 50)


for i in range(0, nentries):

    tree.GetEntry(i)

    muon = []
    antimuon = []
    for lep in range(0, len(tree._lPt)):
        #Check if the lepton indexed lep is a muon. If it is, add index in muon list
        if tree._lpdgId[lep] == 13:
            muon.append(lep)
        #Check if the lepton indexed lep is an antimuon. If it is, add index to antimuon list
        elif tree._lpdgId[lep] == -13:
            antimuon.append(lep)

    if (len(muon) != 2) or (len(antimuon) != 2):
        continue


    dist00 = np.sqrt((tree._lPhi[muon[0]] - tree._lPhi[antimuon[0]]) ** 2
                     + (tree._lEta[muon[0]] - tree._lEta[antimuon[0]]) ** 2)
    dist11 = np.sqrt((tree._lPhi[muon[1]] - tree._lPhi[antimuon[1]]) ** 2
                     + (tree._lEta[muon[1]] - tree._lEta[antimuon[1]]) ** 2)

    dist01 = np.sqrt((tree._lPhi[muon[0]] - tree._lPhi[antimuon[1]]) ** 2
                     + (tree._lEta[muon[0]] - tree._lEta[antimuon[1]]) ** 2)
    dist10 = np.sqrt((tree._lPhi[muon[1]] - tree._lPhi[antimuon[0]]) ** 2
                     + (tree._lEta[muon[1]] - tree._lEta[antimuon[0]]) ** 2)

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
    mass4m = invariantmass4l(tree, 0, 1, 2, 3)

    # rearange pairs so that pair a is always the most massive (~ 90GeV Z boson)
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

    res1 = 10
    res2 = 5

    # S4GeV
    if mass2m_a < 4 or mass2m_b < 4:
        continue

    # ppZZ
    if abs(mass2m_a - 91) < res1 and abs(mass2m_b - 91) < res1:
        continue

    # D20GeV
    #if mass2m_a + mass2m_b < 20:
    #    continue


    pts = []
    pts_max2min = []

    for lep in range(0, len(tree._lPt)):
        pts.append(tree._lPt[lep])
    for lep in range(0, len(tree._lPt)):
        pts_max2min.append(max(pts))
        pts.remove(max(pts))

    #if pts_max2min[2] < 15 or pts_max2min[3] < 10:
    #    continue

    hist_Pt1.Fill(pts_max2min[0])
    hist_Pt2.Fill(pts_max2min[1])
    hist_Pt3.Fill(pts_max2min[2])
    hist_Pt4.Fill(pts_max2min[3])

savehisto(outputfile, hist_Pt1, "hist_Pt_1st", outputdir=outputdirect)
savehisto(outputfile, hist_Pt2, "hist_Pt_2nd", outputdir=outputdirect)
savehisto(outputfile, hist_Pt3, "hist_Pt_3rd", outputdir=outputdirect)
savehisto(outputfile, hist_Pt4, "hist_Pt_4th", outputdir=outputdirect)



