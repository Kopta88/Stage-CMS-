# %%
import os, sys, math, ROOT
from array import array
import numpy as np
from numpy import cos, sin, tan, cosh, sinh
from numpy.linalg import norm
from numpy.ma.core import arccos

from helper import *

root_file = "data.root"
#root_file = "GluGluHToZZTo4L_M125_new.root"

# Pt, Phi, Eta or R
val = "Pt"

# Input file
inputfile = ROOT.TFile(root_file, "READ")

if root_file == "data.root":
    outputdirect = '4MUON/' + val
elif root_file == "GluGluHToZZTo4L_M125_new.root":
    outputdirect = 'Signal/' + val
else:
    print("Invalid root_file, no outputdirect")


# Output file
outputfile = ROOT.TFile("output.root", "RECREATE")

# Accessing the "TTree" from root that contains the data (one entry per event)
tree = inputfile.Get("ntuplizer/tree")

nentries = 5000  # Run only on partial events to speed up the code
nentries = tree.GetEntries()
print("Number of entries: ", nentries)

if val == 'Pt':
    min = 0
    max = 50
elif val == 'R':
    min = 0
    max = 5
else:
    min = -4
    max = 4


hist_val_ma = ROOT.TH1F("hist_"+val+"_ma", ";"+val+" of muon a;Events", 100, min, max)
hist_val_aa = ROOT.TH1F("hist_"+val+"_aa", ";"+val+" of antimuon a;Events", 100, min, max)

hist_val_mb = ROOT.TH1F("hist_"+val+"_ma", ";"+val+" of muon b;Events", 100, min, max)
hist_val_ab = ROOT.TH1F("hist_"+val+"_aa", ";"+val+" of antimuon b;Events", 100, min, max)

hist_val_ma_aa = ROOT.TH2F("hist_"+val+"_pair_a", ";"+val+" of muon;"+val+" of antimuon;Events",
                         50, min, max, 50, min, max)
hist_val_mb_ab = ROOT.TH2F("hist_"+val+"_pair_b", ";"+val+" of muon;"+val+" of antimuon;Events",
                         50, min, max, 50, min, max)

hist_val_Ra = ROOT.TH1F("hist_"+val+"_a", ";"+val+" of pair a;Events", 100, min, max)
hist_val_Rb = ROOT.TH1F("hist_"+val+"_b", ";"+val+" of pair b;Events", 100, min, max)

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

    dist00 = np.sqrt((tree._lPhi[muon[0]]-tree._lPhi[antimuon[0]])**2
                        + (tree._lEta[muon[0]]-tree._lEta[antimuon[0]])**2)
    dist11 = np.sqrt((tree._lPhi[muon[1]]-tree._lPhi[antimuon[1]])**2
                        + (tree._lEta[muon[1]]-tree._lEta[antimuon[1]])**2)

    dist01 = np.sqrt((tree._lPhi[muon[0]]-tree._lPhi[antimuon[1]])**2
                        + (tree._lEta[muon[0]]-tree._lEta[antimuon[1]])**2)
    dist10 = np.sqrt((tree._lPhi[muon[1]]-tree._lPhi[antimuon[0]])**2
                        + (tree._lEta[muon[1]]-tree._lEta[antimuon[0]])**2)

    if (dist00 < dist01) and (dist11 < dist10):
        pair1 = [muon[0], antimuon[0]]
        pair2 = [muon[1], antimuon[1]]
    elif dist00 + dist11 < dist01 + dist10:
        pair1 = [muon[0], antimuon[0]]
        pair2 = [muon[1], antimuon[1]]
    else:
        pair1 = [muon[0], antimuon[1]]
        pair2 = [muon[1], antimuon[0]]

    '''
    if i % 10000 == 0:
        print("\nEvent number", i, "\nNb of reconstructed vertices: ", tree._n_PV)
        print('Pair a')
    for lep in range(0, 2):
        # Every 1000th entry, print some info about the leptons in the event
        if i % 10000 == 0:
            print(
                "Lepton ",
                lep,
                "pt:",
                "%.3f" % tree._lPt[pair1[lep]],
                ", eta:",
                "%.3f" % tree._lEta[pair1[lep]],
                ", phi:",
                "%.3f" % tree._lPhi[pair1[lep]],
                ", pdgid:",
                tree._lpdgId[pair1[lep]],
            )
    if i % 10000 == 0:
        print('Pair b')
    for lep in range(0, 2):
        # Every 1000th entry, print some info about the leptons in the event
        if i % 10000 == 0:
            print(
                "Lepton ",
                lep,
                "pt:",
                "%.3f" % tree._lPt[pair2[lep]],
                ", eta:",
                "%.3f" % tree._lEta[pair2[lep]],
                ", phi:",
                "%.3f" % tree._lPhi[pair2[lep]],
                ", pdgid:",
                tree._lpdgId[pair2[lep]],
                ""
            )
    '''
    mass2m_a = invariantmass2l(tree, pair1[0], pair1[1])
    mass2m_b = invariantmass2l(tree, pair2[0], pair2[1])
    mass4m = invariantmass4l(tree, 0, 1, 2,3)

    # rearange pairs so that pair a is always the most massive (~ 90GeV Z boson)
    if mass2m_a > mass2m_b:
        pair_a = pair1
        pair_b = pair2
    else:
        pair_a = pair2
        pair_b = pair1

    res1 = 10

    # S4GeV
    if mass2m_a < 4 or mass2m_b < 4:
        continue
    # ppZZ
    if abs(mass2m_a-91) < res1 and abs(mass2m_b - 91) < res1:
        continue
    # D20GeV
    #if mass2m_a+mass2m_b < 20 :
    #    continue

    R_a = np.sqrt((tree._lPhi[pair_a[0]]-tree._lPhi[pair_a[1]])**2
                        + (tree._lEta[pair_a[0]]-tree._lEta[pair_a[1]])**2)
    R_b = np.sqrt((tree._lPhi[pair_b[0]]-tree._lPhi[pair_b[1]])**2
                        + (tree._lEta[pair_b[0]]-tree._lEta[pair_b[1]])**2)

    if val == 'Pt':
        hist_val_ma.Fill(tree._lPt[pair_a[0]])
        hist_val_aa.Fill(tree._lPt[pair_a[1]])
        hist_val_mb.Fill(tree._lPt[pair_b[0]])
        hist_val_ab.Fill(tree._lPt[pair_b[1]])

        hist_val_ma_aa.Fill(tree._lPt[pair_a[0]], tree._lPt[pair_a[1]])
        hist_val_mb_ab.Fill(tree._lPt[pair_b[0]], tree._lPt[pair_b[1]])
    elif val == 'Phi':
        hist_val_ma.Fill(tree._lPhi[pair_a[0]])
        hist_val_aa.Fill(tree._lPhi[pair_a[1]])
        hist_val_mb.Fill(tree._lPhi[pair_b[0]])
        hist_val_ab.Fill(tree._lPhi[pair_b[1]])

        hist_val_ma_aa.Fill(tree._lPhi[pair_a[0]], tree._lPhi[pair_a[1]])
        hist_val_mb_ab.Fill(tree._lPhi[pair_b[0]], tree._lPhi[pair_b[1]])
    elif val == 'Eta':
        hist_val_ma.Fill(tree._lEta[pair_a[0]])
        hist_val_aa.Fill(tree._lEta[pair_a[1]])
        hist_val_mb.Fill(tree._lEta[pair_b[0]])
        hist_val_ab.Fill(tree._lEta[pair_b[1]])

        hist_val_ma_aa.Fill(tree._lEta[pair_a[0]], tree._lEta[pair_a[1]])
        hist_val_mb_ab.Fill(tree._lEta[pair_b[0]], tree._lEta[pair_b[1]])
    elif val == 'R':
        hist_val_Ra.Fill(R_a)
        hist_val_Rb.Fill(R_b)
    else:
        print("Incorrect string value for 'val' ")

savehisto(outputfile, hist_val_ma, "hist_"+val+"_Am", outputdir=outputdirect)
savehisto(outputfile, hist_val_aa, "hist_"+val+"_Aa", outputdir=outputdirect)
savehisto(outputfile, hist_val_mb, "hist_"+val+"_Bm", outputdir=outputdirect)
savehisto(outputfile, hist_val_ab, "hist_"+val+"_Ba", outputdir=outputdirect)

#savehisto(outputfile, hist_val_ma_aa, "hist_"+val+"_m_a_A", outputdir=outputdirect)
#savehisto(outputfile, hist_val_mb_ab, "hist_"+val+"_m_a_B", outputdir=outputdirect)

#savehisto(outputfile, hist_val_Ra, "hist_"+val+"_a", outputdir=outputdirect)
#savehisto(outputfile, hist_val_Rb, "hist_"+val+"_b", outputdir=outputdirect)





