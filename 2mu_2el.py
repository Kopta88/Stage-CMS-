# %%
import os, sys, math, ROOT
from array import array
import numpy as np
from numpy import cos, sin, tan, cosh, sinh
from numpy.linalg import norm
from numpy.ma.core import arccos

from helper import *

# Input file
inputfile = ROOT.TFile("data.root", "READ")
outputdirect = '2Muon2Elect'
# Output file
outputfile = ROOT.TFile("output.root", "RECREATE")

# Accessing the "TTree" from root that contains the data (one entry per event)
tree = inputfile.Get("ntuplizer/tree")

nentries = 1000  # Run only on partial events to speed up the code
nentries = tree.GetEntries()
print("Number of entries: ", nentries)

hist_mass4l = ROOT.TH1F("hist_mass4leptons", ";Invariant Mass of four leptons;Events", 100, 0, 20)
hist_mass2m = ROOT.TH1F("hist_massmuonpair", ";Invariant Mass of muon pair;Events", 100, 0, 200)
hist_mass2e = ROOT.TH1F("hist_masselectronpair", ";Invariant Mass of electron pair;Events", 100, 0, 200)
hist_mass2m_2e = ROOT.TH2F("hist_massleptonpairs", ";Invariant Mass of muon pair;Invariant Mass of electron pair;Events",
                         100, 0, 200, 100, 0, 200)
hist_mass4l_2m = ROOT.TH2F("hist_massmuonpair1&mass4l", ";4 lepton Invariant Mass;Invariant Mass of muon pair;Events",
                         100, 0, 500, 100, 0, 200)
hist_mass4l_2e = ROOT.TH2F("hist_masselectronpair2&mass4l", ";4 lepton Invariant Mass;Invariant Mass of electron pair;Events",
                         100, 0, 500, 100, 0, 200)


for i in range(0, nentries):

    tree.GetEntry(i)

    pdgId_absum = 0
    pdgId_sum = 0
    for lep in range(0, len(tree._lPt)):
        pdgId_sum += tree._lpdgId[lep]
        pdgId_absum += abs(tree._lpdgId[lep])

    # first selector skips the entry if not 2 muon-antimuon pairs
    if pdgId_sum != 0 or pdgId_absum != 48:
        continue

    muons = []
    electrons = []
    for lep in range(0, len(tree._lPt)):
        if abs(tree._lpdgId[lep]) == 13:
            muons.append(lep)
        else:
            electrons.append(lep)

    mass2mu = invariantmass2l(tree, muons[0], muons[1])
    mass2el = invariantmass2l(tree, electrons[0], electrons[1])
    mass4l = invariantmass4l(tree, 0, 1, 2,3)

    #if mass2mu < 4 or mass2el < 4:
    #    continue

    # print("invariant mass of first pair: ", mass2l_a)
    # print("invariant mass of second pair: ", mass2l_b, '\n')

    hist_mass4l.Fill(mass4l)
    #hist_mass2m.Fill(mass2mu)
    #hist_mass2e.Fill(mass2el)
    #hist_mass2m_2e.Fill(mass2mu, mass2el)
    #hist_mass4l_2m.Fill(mass4l, mass2mu)
    #hist_mass4l_2e.Fill(mass4l, mass2el)


#savehisto(outputfile, hist_mass4l, "hist_mass4l_x20", outputdir=outputdirect)
#savehisto(outputfile, hist_mass4l_2m, "hist_mass4l_2m", outputdir=outputdirect)
#savehisto(outputfile, hist_mass4l_2e, "hist_mass4l_2e_no4GeV", outputdir=outputdirect)
#savehisto(outputfile, hist_mass2m_2e, "hist_mass2m_2e_no4GeV", outputdir=outputdirect)
savehisto(outputfile, hist_mass2m, "hist_mass2m", outputdir=outputdirect)
savehisto(outputfile, hist_mass2e, "hist_mass2e", outputdir=outputdirect)




