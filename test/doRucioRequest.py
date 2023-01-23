import subprocess
import sys
import os

def bash( bashCommand ):
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	#process = subprocess.Popen(bashCommand.split())
	output, error = process.communicate()
	return output ,error

def bashout( command ):
	output = subprocess.check_output( command, shell=True)
	return output	

def doCommand( command ):
	output = os.system( command )
	return output

gmsb_list = [
#    'cms:/GMSB_L-100TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    'cms:/GMSB_L-400TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-0p1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-0p1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    'cms:/GMSB_L-600TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
]

gjets_list = [
    'cms:/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-4cores5k_106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    'cms:/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    'cms:/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    'cms:/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    'cms:/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
]

#runlist = gjets_list
runlist = gmsb_list
for thefile in runlist :

	thecommand = 'rucio add-rule --ask-approval --lifetime 5184000 ' + thefile + ' 1 T3_US_FNALLPC'
	print( '---------------------------------------------------' )
	print( thecommand )
	doCommand( thecommand )
	print( '---------------------------------------------------' )

#print( bashout( 'hadd llpgana_HTo2LongLivedTo4b_t37MC_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root ' + theFileList ) )	
