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
options.register('globalTag','106X_dataRun2_v28',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
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
outfilename = 'llpgana_mc_AODSIM_GMSB17_ST.root'

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
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 10000

## Define the input source
aodpath_1k_450_100k = '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_14TeV_pythia8/AODSIM/110X_mcRun3_2021_realistic_v6-v2/'
aodpath_1k_450_10k = '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV_pythia8/AODSIM/110X_mcRun3_2021_realistic_v6-v2/'
aodpath_125_25_15k = '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-15000mm_TuneCP5_14TeV_pythia8/AODSIM/110X_mcRun3_2021_realistic_v6-v2/'
aodpath21_1k_450_100k = '/store/mc/Run3Summer21DRPremix/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV-pythia8/AODSIM/120X_mcRun3_2021_realistic_v6-v2/'

lpcpath_350_600 = 'file:/eos/uscms/store/mc/RunIIFall17DRPremix/GMSB_L-350TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

		  #'file:jwk_reco_data_DIGI2RAW.root'),

         ## HTo2LongLivedTo4b

        #'/store/mc/Run3Winter21DRMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_14TeV-pythia8/MINIAODSIM/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280000/17bd2d71-8a76-46c5-947a-7ea2b1df44b6.root'

		# AOD

		 #'file:/eos/uscms/store/mc/RunIIFall17DRPremix/GMSB_L-100TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/2C55A98D-E4D7-E811-AC40-002590491B1E.root'
		 lpcpath_350_600+'120000/80762156-99D6-E811-8942-34E6D7E3879B.root'

         #'file:BCB550D6-CAB3-5C4F-8866-77897305A646.root',#aodpath_125_25_15k
         #cmssw12XX only #'file:bc04e7b9-31c7-4bec-a396-258ba40b8bd5.root',#aodpath21_1k_450_100k

         #aodpath_1k_450_100k+'10000/07C1D360-FDE3-B04C-8E7A-DDE4275C7F04.root',
         #aodpath_1k_450_100k+'10000/2EBC0785-9656-2B43-A9D2-FFEB5032E66A.root',
         #aodpath_1k_450_100k+'10000/4002AC07-8CA1-0643-9237-F38929578E9D.root',
         #aodpath_1k_450_100k+'10000/447E972B-3F0E-E646-AB2E-8D1D530CFE88.root',
         #aodpath_1k_450_100k+'10000/6C3EEB7D-B4DD-2E4C-AA3C-8BFA5109627E.root',
         #aodpath_1k_450_100k+'10000/701683AD-0961-B346-98C4-349653C33F5C.root',
         #aodpath_1k_450_100k+'10000/7080FB51-EF47-D74F-B68A-187AFCF0278A.root',
         #aodpath_1k_450_100k+'10000/7BE0C49A-88F0-0E4E-9CB5-06DD4CBEB538.root',
         #aodpath_1k_450_100k+'10000/DF28A0DC-1D69-2147-90FB-546BD0443684.root',
         #aodpath_1k_450_100k+'10000/FFC58857-0181-8D41-BB92-1DBA57796135.root',
         #aodpath_1k_450_100k+'240000/49A0EB69-E794-2D4D-9955-4A3ECC27B6F3.root',

         #aodpath_1k_450_10k+'10000/06A10966-0ABE-A94D-A5C0-BA43FB64B0AB.root',
         #aodpath_1k_450_10k+'10000/190B92B5-4416-A346-B0A9-B9980F15563C.root',
         #aodpath_1k_450_10k+'10000/40AF3AF3-6F15-1341-BF57-6A7C2ACBD01F.root',
         #aodpath_1k_450_10k+'10000/4DA81C18-5BCD-5E44-B7AA-AE4963A8D1A5.root',
         #aodpath_1k_450_10k+'10000/61B67C46-6181-EA42-90B6-6680A697CF1B.root',
         #aodpath_1k_450_10k+'10000/7FEE9823-303D-0D4C-A852-97CBA777CBFA.root',
         #aodpath_1k_450_10k+'10000/A449E283-46F7-E74D-B6BF-15862B85F3D8.root',
         #aodpath_1k_450_10k+'10000/D0488318-FB3D-7041-BB72-A009D26248E5.root',
         #aodpath_1k_450_10k+'10000/F9C14E6C-6FC0-F14A-BA16-6E44B7480009.root',
         #aodpath_1k_450_10k+'50000/4E0AC838-A6FF-B549-8441-A86B04B2BEFA.root',

         #aodpath_125_25_15k+'240000/2543A8DF-4540-FE45-A496-9EF8D7918E35.root',
         #aodpath_125_25_15k+'240000/6299C20E-A359-8B4F-A4C3-858E50BC461E.root',
         #aodpath_125_25_15k+'240000/6EA51191-B27D-B547-920E-66EF94292870.root',
         #aodpath_125_25_15k+'240000/74F70EBA-B107-9348-B707-5FA925F5B13F.root',
         #aodpath_125_25_15k+'240000/8E42896A-57F9-284A-A1A3-26857DD2A4F5.root',
         #aodpath_125_25_15k+'240000/9381F9E2-2338-5049-AB04-0EDDD054F299.root',
         #aodpath_125_25_15k+'240000/A330CAF2-7A94-4F40-B4DE-2287AB3E5831.root',
         #aodpath_125_25_15k+'240000/B4282E0A-8B1A-D94C-84CB-A85F0F210D1E.root',
         #aodpath_125_25_15k+'240000/B5CECD56-5DA4-E44C-B5D2-508043F27FCE.root',
         #aodpath_125_25_15k+'240000/BCB550D6-CAB3-5C4F-8866-77897305A646.root',
         #aodpath_125_25_15k+'240000/C7D370A7-A4B4-A442-BB53-4A5F0EDF10AE.root',

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
		
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/059EA86A-6D43-074E-B511-1B33CFBBA6F3.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/06891533-FC18-7646-AAEF-5707CD0C303C.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/06DF5368-B816-114F-B363-CBB4D6A32E28.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0715EFFE-09C9-7347-B37E-CE5CCE65B80A.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0738A4E6-52AC-3142-9CB7-4E059E55D4B4.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/074128AA-4BEC-6C45-9652-F6C93CF64D8E.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/08686EAC-2BF0-9848-8161-F926A417C735.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/087313BD-CAC2-DA4B-B0E8-65B3EC9C100D.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0886CEBB-2006-3140-AFE1-722ACD2E79F1.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/08B336D5-3BD6-CB43-9B87-311488546C76.root',

        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/09181BBE-2688-BF45-9210-E5272F847646.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/094109BF-B509-2E41-8962-A828AD2956B1.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/09413E06-E33D-4941-BB6A-B38E2D3710B0.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/098AF77C-D5BB-4348-AAFB-E186E1730D28.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/09C67C38-22ED-BB44-A84C-DA1E8AC42EE2.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0A0F0808-5DA0-474F-928A-30D916D9AB1B.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0A529A9C-C6DC-B04E-BC81-A64144461C4F.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0A75B49C-CA2D-B249-8226-0EBE36E6B7F9.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0BA377C2-8CA6-254D-889B-789E6AF31377.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0C9F9B98-3AE4-4F46-A1B6-59BAAC2A0A34.root',

        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0CB669FF-60C0-F84F-8165-6D75030DE776.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0CC3059D-4737-B349-8C93-F8330263F87E.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0D0AF398-45F9-4C47-B639-DE12528362BC.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0D2A5724-E01E-CC4B-A76F-1152B4AE56D9.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0DB2F6E9-9E2F-0B48-BC9A-75620901E60D.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0E7EC6AF-3277-CE45-A609-B7186BC7FB7C.root',

        ),##<<>>fileNames = cms.untracked.vstring
)##<<>>process.source = cms.Source("PoolSource",


## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))#ST
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))#T
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))#LT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500))#US
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(12500))#VS
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(25000))#S
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))#SM
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250000))#M
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500000))#L
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))#F

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
   ## electrons
   #electrons = cms.InputTag("slimmedElectrons"),
   electrons = cms.InputTag("gedGsfElectrons"),
   ## muons
   #muons = cms.InputTag("slimmedMuons"),
   muons = cms.InputTag("muons"),
   ## photons
   #gedPhotons = cms.InputTag("slimmedPhotons"),
   gedPhotons = cms.InputTag("gedPhotons"),
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

