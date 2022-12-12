import os, re
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

## Flags
options.register('hasGenInfo',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to get pcalo in mc');

## object prep cuts
#options.register('jetpTmin',15.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'jet pT minimum cut');
#options.register('jetEtamax',3.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'jet eta maximum cut');
#options.register('jetIDmin',1,VarParsing.multiplicity.singleton,VarParsing.varType.int,'jet ID minimum cut');
#options.register('rhEmin',1.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'recHit energy minimum cut');
#options.register('phpTmin',20.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'photon pT minimum cut');
#options.register('phIDmin','none',VarParsing.multiplicity.singleton,VarParsing.varType.string,'photon ID minimum cut');

## lepton prep cuts
#options.register('ellowpTmin',20.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'electron low pT min cut');
#options.register('elhighpTmin',50.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'electron high pT min cut');
#options.register('mulowpTmin',20.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'muon low pT minimum cut');
#options.register('muhighpTmin',50.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'muon high pT minimum cut');

## GT to be used
options.register('globalTag','94X_mc2017_realistic_v14',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','106X_dataRun2_v28',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','112X_mcRun3_2021_realistic_v16',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#112X_mcRun3_2021_realistic_v16

## processName
options.register('processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,'process name to be considered');

## outputFile Name
#outfilename = 'llpgana_t35a_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root'
#outfilename = 'llpgana_t36L_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 35a + ph time
#outfilename = 'llpgana_t36L_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 35a + ph time
#outfilename = 'llpgana_t37MC_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 35a + ph time
#outfilename = 'llpgana_t68MC_eigen_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # adding genJet Parton info
#outfilename = 'llpgana_mc_t75M_pheigen95t60r9_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 73 w/ 3DProfile & sphcalrot : ebp ebn + ltime sum flip
#outfilename = 'llpgana_mc_t76SM_pheigen95t60r9_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as May 16 2022 update with energy/delaytime info
#outfilename = 'llpgana_mc_t79M_pheigen95t60r9_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 76 + ootPhotons only
#outfilename = 'llpgana_mc_t80M_pheigen95t60r9_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 79 + genJetDrmatch
#outfilename = 'llpgana_mc_AODSIM_GMSB100TeV_v80a_selected_FL.root'
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v3.root' # finishd addig rechits
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v4.root' # testing saved rhIDsGroup for photons
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v5.root' # testing valuemap id cut loose for photons
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v6.root' # adding photon information
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v7.root' # adding ootphoton information
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v8.root' #  added genparticle info for photons, changed genjet to best dr match
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v10.root' # modded genpart to match only photons, mom count on genpart collection
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v11.root' # combined oot & ged phtons into single collection, modded gen particle llp is coding
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v12.root' # modifed genpart to look for susy llp particles, included met vars in output
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v13.root' # added genprt collection with phogen reffrenced to it
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v14.root' # added corrected met and inserted jet SC rhID groups
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v15.root' # corrected eta and phi for oot photon gen match
#outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v16.root' # eta/phi correction test for oot phos
outfilename = 'llpgana_mc_AODSIM_ntuplizer_test_v17.root' # turned of jets, ect., changed genpart stablity check, slimmed output

options.register('outputFileName',outfilename,VarParsing.multiplicity.singleton,VarParsing.varType.string,'output file name created by cmsRun');

## parsing command line arguments
options.parseArguments()

## Define the CMSSW process
from Configuration.StandardSequences.Eras import eras
process = cms.Process(options.processName,eras.Run2_2018)

## Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

#process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
#process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Configuration.StandardSequences.EndOfProcess_cff')

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 10000

## Define the input source
aodpath_1k_450_100k = '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_14TeV_pythia8/AODSIM/110X_mcRun3_2021_realistic_v6-v2/'
aodpath_1k_450_10k = '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV_pythia8/AODSIM/110X_mcRun3_2021_realistic_v6-v2/'
aodpath_125_25_15k = '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-15000mm_TuneCP5_14TeV_pythia8/AODSIM/110X_mcRun3_2021_realistic_v6-v2/'
aodpath21_1k_450_100k = '/store/mc/Run3Summer21DRPremix/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV-pythia8/AODSIM/120X_mcRun3_2021_realistic_v6-v2/'

lpcpath_350_600 = 'file:/eos/uscms/store/mc/RunIIFall17DRPremix/GMSB_L-350TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/'
lpcpath_GMSB = 'file:/eos/uscms/store/mc/RunIIFall17DRPremix/'
gmsbaodsim = '_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/'
gmsbaodsim2 = '_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

		'root://cmsxrootd-site.fnal.gov//store/mc/RunIIFall17DRPremix/GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/80000/58DCC006-3FB2-E811-94AE-AC1F6B0DE454.root',	

		  #'file:jwk_reco_data_DIGI2RAW.root'),

         ## HTo2LongLivedTo4b

        #'/store/mc/Run3Winter21DRMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_14TeV-pythia8/MINIAODSIM/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280000/17bd2d71-8a76-46c5-947a-7ea2b1df44b6.root'

		# AOD

		 #'file:/eos/uscms/store/mc/RunIIFall17DRPremix/GMSB_L-100TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/2C55A98D-E4D7-E811-AC40-002590491B1E.root'
		 #lpcpath_350_600+'120000/80762156-99D6-E811-8942-34E6D7E3879B.root',
         #lpcpath_350_600+'120000/322875DC-DDD6-E811-8C5F-001E675A68C4.root',
         #lpcpath_350_600+'120000/603C58F0-A9D6-E811-8907-E0071B74AC00.root',
         #lpcpath_350_600+'120000/66DB15F0-DDD6-E811-A7A3-90E2BACBAD58.root',
         #lpcpath_350_600+'120000/8A9DEBDE-DDD6-E811-838E-D4AE526DF2E1.root',
         #lpcpath_350_600+'120000/A8DE3FA7-99D6-E811-92C5-34E6D7BEAF0E.root'

		 # AODSIM GMSB model		

         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10000cm'+gmsbaodsim2+'10000/162DBEEE-DC29-E911-843B-0CC47A745294.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10000cm'+gmsbaodsim2+'10000/BCE6A6F2-A929-E911-AEA9-24BE05C63651.root',
         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10000cm'+gmsbaodsim2+'100000/22F81F5B-3D33-E911-906C-0CC47AD24D28.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10000cm'+gmsbaodsim2+'100000/3CDA2489-9132-E911-B4E3-0025905D1E08.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10000cm'+gmsbaodsim2+'100000/44349E62-F131-E911-97A6-848F69FD09D7.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1000cm'+gmsbaodsim+'270000/26A7404E-C5DA-E811-92D3-001E675A6AB8.root',
         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1000cm'+gmsbaodsim+'270000/66971865-6ADA-E811-A9C2-002590D9D8C0.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1000cm'+gmsbaodsim+'270000/6E716649-6ADA-E811-B9E5-0025901AC0FC.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1000cm'+gmsbaodsim+'270000/A89ECEA5-C5DA-E811-AB4B-D4AE526DDB3F.root',
         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10cm'+gmsbaodsim+'80000/08A3D920-E5B2-E811-85C5-0CC47A4DEEE4.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10cm'+gmsbaodsim+'80000/0C51F68D-E5B2-E811-8703-002590E39F2E.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10cm'+gmsbaodsim+'80000/2097D13F-30B3-E811-9FB3-0025905C96EA.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1200cm'+gmsbaodsim+'60000/605B95ED-8AD8-E811-B640-001E67A3F70E.root',
         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1200cm'+gmsbaodsim+'60000/68C5DD84-EAD6-E811-9028-44A842CFD60C.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1200cm'+gmsbaodsim+'60000/846C6D28-D4D8-E811-A8BE-00266CFFBEB4.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1200cm'+gmsbaodsim+'60000/8EC1497F-E9D6-E811-A390-44A842CFC98B.root',
         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-200cm'+gmsbaodsim+'120000/ACDBAA3E-9ED8-E811-813B-1866DA89095D.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-200cm'+gmsbaodsim+'120000/AEA9923B-9ED8-E811-87D6-34E6D7BDDECE.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-400cm'+gmsbaodsim+'80000/301550E4-81B5-E811-8362-FA163EF08F5B.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-400cm'+gmsbaodsim+'80000/34D601E8-3CB2-E811-B2BA-24BE05C6B701.root',
         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-400cm'+gmsbaodsim+'80000/4ACDB9DD-C0AF-E811-A364-24BE05C6D731.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-400cm'+gmsbaodsim+'80000/7411AE64-81B5-E811-9F51-FA163E301A96.root',
         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-600cm'+gmsbaodsim+'60000/E25EAAE9-A5D6-E811-8502-34E6D7E05F0E.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-600cm'+gmsbaodsim+'60000/F65FFA8E-09D7-E811-81E5-A4BF0102A4F5.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-800cm'+gmsbaodsim+'270000/B624F4E1-86D8-E811-BD73-A0369F83633E.root',
         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-800cm'+gmsbaodsim+'270000/DE490468-30D8-E811-99AC-0025B3E015D2.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-800cm'+gmsbaodsim+'80000/2C3F963F-3BB2-E811-BA81-D8D385AF8902.root',

         #'file:BCB550D6-CAB3-5C4F-8866-77897305A646.root',#aodpath_125_25_15k
         #cmssw12XX only #'file:bc04e7b9-31c7-4bec-a396-258ba40b8bd5.root',#aodpath21_1k_450_100k

		 # HTo2LongLivedTo4b

         #aodpath_1k_450_100k+'10000/07C1D360-FDE3-B04C-8E7A-DDE4275C7F04.root',
         #aodpath_1k_450_100k+'10000/2EBC0785-9656-2B43-A9D2-FFEB5032E66A.root',
         #aodpath_1k_450_100k+'10000/4002AC07-8CA1-0643-9237-F38929578E9D.root',

         #aodpath_1k_450_10k+'10000/06A10966-0ABE-A94D-A5C0-BA43FB64B0AB.root',
         #aodpath_1k_450_10k+'10000/190B92B5-4416-A346-B0A9-B9980F15563C.root',
         #aodpath_1k_450_10k+'50000/4E0AC838-A6FF-B549-8441-A86B04B2BEFA.root',

         #aodpath_125_25_15k+'240000/2543A8DF-4540-FE45-A496-9EF8D7918E35.root',
         #aodpath_125_25_15k+'240000/6299C20E-A359-8B4F-A4C3-858E50BC461E.root',
         #aodpath_125_25_15k+'240000/6EA51191-B27D-B547-920E-66EF94292870.root',

		  ## EGamma

        #'/store/data/Run2018D/EGamma/MINIAOD/22Jan2019-v2/70001/F0FE59FC-F29E-904B-A1BC-817C1CB09A7E.root'

        #'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/062A29F3-416F-6A43-8E0A-90BE80607677.root',
        #'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/0632DF20-D859-EE4B-8FC4-3ECC56B2987D.root',
        #'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/0685BE81-54CF-3A4C-8D66-8E0E059D263F.root',
        #'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/06CBEAC6-F3C1-BE48-96D9-A7EF3CF42048.root',
        #'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/06CC5637-77D3-8145-B3BD-0AF5F16410AB.root',
        #'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/06DE38F0-4330-C546-8D1A-0865C50FAA14.root',
        #'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/071F2728-A0CA-9446-A1E3-38B4FB08D6AE.root',
        #'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/077E7C54-C588-F642-BBA8-0960D5FD7A96.root',
        #'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/07971C26-63F0-D943-BA77-88D6F636F780.root',
        #'/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/07A40762-F4E1-9E4D-8750-7DD65EA01B8A.root',

		  ## JetHT

        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/00B87525-94D1-C741-9B03-00528106D15A.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/00F0369F-B9FA-5643-B0F4-F286B80B30B5.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0173DB9C-3E46-4741-B80B-9671E9E917CA.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/017C1808-933C-0D48-95B5-A8B9A0B16601.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/01C71814-B84D-1D4A-8696-61E404084059.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/020C4D59-9FF4-BA4C-B24C-6B77F9F0CC20.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/027298C8-F7C2-7848-B241-DD6F07F41C2F.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0388A30D-E3A0-2947-8F6A-52B8B3F4C381.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/04939F8B-6606-3F4C-92BF-6D7D950376DC.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/054042A4-058C-194E-9756-038004285B82.root',
		
        ),##<<>>fileNames = cms.untracked.vstring
)##<<>>process.source = cms.Source("PoolSource",


## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))#ST
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))#TTi
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))#LT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500))#US
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(12500))#VS
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(25000))#SM
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))#MS
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250000))#MD
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500000))#LG
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))#FL

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = options.globalTag  

## Create output file
## Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", 
		                   fileName = cms.string(options.outputFileName))
				   
# Make the tree 
process.tree = cms.EDAnalyzer("LLPgammaAnalyzer_AOD",
   ## flags
   hasGenInfo = cms.bool(options.hasGenInfo),
   ## additional collections
   ## tracks
   #tracks = cms.InputTag("unpackedTracksAndVertices"),
   tracks = cms.InputTag("generalTracks"),
   ## vertices
   #vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   vertices = cms.InputTag("offlinePrimaryVertices"),
   ## pfcandidates
   #pfcandidates = cms.InputTag("packedPFCandidates"),
   pfcandidates = cms.InputTag("particleFlow"),
   particleflow = cms.InputTag("particleFlow",""),	
   pfcanphomap = cms.InputTag("particleFlow","photons"),
   pfcanootphomap = cms.InputTag("particleFlow","photons"),
   pfcanelemap = cms.InputTag("particleFlow","electrons"),
   ## rho
   #rho = cms.InputTag("fixedGridRhoFastjetAll"), #fixedGridRhoAll
   rho = cms.InputTag("fixedGridRhoAll"),
   ## METs
   #mets = cms.InputTag("slimmedMETs"),
   mets = cms.InputTag("pfMet"),
   ## jets
   #jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
   #jets = cms.InputTag("slimmedJets"),
   jets = cms.InputTag("ak4PFJets"),
   genjets = cms.InputTag("ak4GenJets",""),
   calojets = cms.InputTag("ak4CaloJets",""),
   ## electrons
   #electrons = cms.InputTag("slimmedElectrons"),
   electrons = cms.InputTag("gedGsfElectrons"),
   ## muons
   #muons = cms.InputTag("slimmedMuons"),
   muons = cms.InputTag("muons"),
   ## photons
   #gedPhotons = cms.InputTag("slimmedPhotons"),
   #gedPhotons = cms.InputTag("gedPhotons"),
   gedPhotons = cms.InputTag("photons"),
   phoCBIDLooseMap = cms.InputTag("PhotonIDProd", "PhotonCutBasedIDLooseEM"), 
   #ootPhotons = cms.InputTag("slimmedOOTPhotons"),
   ootPhotons = cms.InputTag("ootPhotons"),
   ## ecal recHits
   #recHitsEB = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
   #recHitsEE = cms.InputTag("reducedEgamma", "reducedEERecHits"),
   recHitsEB = cms.InputTag("reducedEcalRecHitsEB"),
   recHitsEE = cms.InputTag("reducedEcalRecHitsEE"),
   ## superclusters
   #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters"),
   superClusters = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALBarrel"),
   #ootSuperClusters = cms.InputTag("reducedEgamma", "reducedOOTSuperClusters"),
   ootSuperClusters = cms.InputTag("particleFlowSuperClusterOOTECAL", "particleFlowSuperClusterOOTECALBarrel"), 
   ## caloclusters
   #caloClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters"),
   caloClusters = cms.InputTag("particleFlowEGamma", "EBEEClusters"),
   ## gen info
   genEvt = cms.InputTag("generator", ""),
   gent0 = cms.InputTag("genParticles", "t0"), 
   genxyz0 = cms.InputTag("genParticles", "xyz0"), 
   pileups = cms.InputTag("addPileupInfo", ""),
   #pileups = cms.InputTag("mixData", ""),
   genParticles = cms.InputTag("genParticles", ""),		

)##<<>>process.tree = cms.EDAnalyzer("LLPgammaAnalyzer_aod"


# Set up the path
process.tree_step = cms.EndPath(process.tree)
process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.content_step = cms.Path(process.content)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.schedule = cms.Schedule(process.tree_step)
process.options = cms.untracked.PSet()

#do not add changes to your config after this point (unless you know what you are doing)
#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#process=convertToUnscheduled(process)

# customisation of the process.
#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
#process = miniAOD_customizeAllData(process)
# End of customisation functions

process.options = cms.untracked.PSet( 
    #numberOfThreads = cms.untracked.uint32(4), 
    #numberOfStreams = cms.untracked.uint32(4), 
    SkipEvent = cms.untracked.vstring('ProductNotFound'), 
    #wantSummary = cms.untracked.bool(True) 
)

