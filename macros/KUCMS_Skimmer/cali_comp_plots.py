#  jack w king 3

#from numpy import array
from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

from jwk_tdr_style_py import *

def dostack():

	fOutFile = TFile( "cali_QCD_Gjets_phorhE5_ewmo_diff_hists.root", "RECREATE" )

	diffhist = TH2F("DiffMap","DiffMap",171,-85.5,85.5,360,0.5,360.5)
	disthist = TH1F("DiffDist","DiffDist",240,-0.6,0.6) 
	etadiffhist = TH1F("EtaDiffHist","EtaDiffHist",171,-85.5,85.5)
	etadist = TH1F("EtaDiffDist","EtaDiffDist",240,-0.06,0.06)

	tree = ""
	histname = "AveXtalRatioRecTimeEBMap"
	etahist = "AveEtaRecTimeEta"
	infile1 = "cali_root_files/KUCMS_QCD_R17_v16_phorhE5_ewmo_Cali.root"
	infile2 = "cali_root_files/KUCMS_GJets_R17_v16_phorhE5_ewmo_Cali.root"    

	f1 = TFile.Open(infile1)
	if tree == '' : hist = histname
	else : hist = tree+'/'+histname

	f2 = TFile.Open(infile2)
	if tree == '' : hist = histname
	else : hist = tree+'/'+histname

	orighist1 = f1.Get(hist)
	orighist2 = f2.Get(hist)
	etahist1 = f1.Get(etahist)
	etahist2 = f2.Get(etahist)

	for eta in range(1,171):

		ebinval1 = float(etahist1.GetBinContent(eta))
		ebinval2 = float(etahist2.GetBinContent(eta))
		ediff = ebinval1 - ebinval2
		ebinerr = sqrt(float(etahist1.GetBinError(eta))*float(etahist1.GetBinError(eta)) + float(etahist2.GetBinError(eta))*float(etahist2.GetBinError(eta)))
		etadiffhist.SetBinContent(eta,ediff)
		etadiffhist.SetBinError(eta,ebinerr);
		etadist.Fill(ediff)

		for phi in range(1,360):

			binval1 = float(orighist1.GetBinContent(eta,phi))
			binval2 = float(orighist2.GetBinContent(eta,phi))
			diffbin = binval1 - binval2
			diffhist.SetBinContent(eta,phi,diffbin)
			disthist.Fill(diffbin)

	fOutFile.cd()
	diffhist.Write()
	disthist.Write()
	etadiffhist.Write()
	etadist.Write()

dostack()

