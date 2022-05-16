import os, re
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

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

## processName
options.register('processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,'process name to be considered');

## outputFile Name
#outfilename = 'llpgana_t35a_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root'
#outfilename = 'llpgana_t36L_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 35a + ph time
#outfilename = 'llpgana_t36L_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 35a + ph time
#outfilename = 'llpgana_t39L_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 35a + ph time
#outfilename = 'llpgana_t47S_eigen_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 35a + ph time
#outfilename = 'llpgana_t55S_eigen_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 47 with eigins ( spherical and ieipt in units of time )
#outfilename = 'llpgana_t56S_eigen_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 56+ eigns with wted spherical ( reso*dt*dt + angle+PI if dt negitive )
#outfilename = 'llpgana_t68L_eigen_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 56 with res wted 2d histo
#outfilename = 'llpgana_t69L_eigen_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 68 with ltsum for spher and oval egin in ieipt
#outfilename = 'llpgana_t72L_eigen_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 69 with code if/then loop bug corrected
#outfilename = 'llpgana_t73S_eigen95t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 72 w/ sphcty cut &/or #rh cut for sphcal eta/time profl
#outfilename = 'llpgana_t74L_eigen95t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 73 w/ 3DProfile & ignore 1st couple of etas for sphcal
#outfilename = 'llpgana_t75L_pheigen95t60r9_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 73 w/ 3DProfile & sphcal rot :  ebp ebn + ltime sum for flip
#outfilename = 'llpgana_t74M_pheigen95t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root'#
#outfilename = 'llpgana_t75M_eigen95t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # added 2sigma filter for eta time profile (exclude t>1 eta t map)
#outfilename = 'llpgana_t79L_eigen95t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # add by cluster slope extraction
#outfilename = 'llpgana_t81M_eigen95t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # add by cluster slope extraction
#outfilename = 'llpgana_t82M_eigen95t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 81 + crystalball
#outfilename = 'llpgana_t84L_eigen95t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 81 + 2sided crystalball + slope algiment for 2D
#outfilename = 'llpgana_t85M_eigen95t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 81 + 2sided crystalball + slope algiment for 2D + slope cuts
outfilename = 'llpgana_t88M_eigen96t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 81 + 2sided crystalball + slope algiment for 2D + slope cuts

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

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Configuration.StandardSequences.EndOfProcess_cff')

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

## Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

		  #'file:jwk_reco_data_DIGI2RAW.root'),

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

        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/00B87525-94D1-C741-9B03-00528106D15A.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/00F0369F-B9FA-5643-B0F4-F286B80B30B5.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0173DB9C-3E46-4741-B80B-9671E9E917CA.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/017C1808-933C-0D48-95B5-A8B9A0B16601.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/01C71814-B84D-1D4A-8696-61E404084059.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/020C4D59-9FF4-BA4C-B24C-6B77F9F0CC20.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/027298C8-F7C2-7848-B241-DD6F07F41C2F.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0388A30D-E3A0-2947-8F6A-52B8B3F4C381.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/04939F8B-6606-3F4C-92BF-6D7D950376DC.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/054042A4-058C-194E-9756-038004285B82.root',
		
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/059EA86A-6D43-074E-B511-1B33CFBBA6F3.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/06891533-FC18-7646-AAEF-5707CD0C303C.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/06DF5368-B816-114F-B363-CBB4D6A32E28.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0715EFFE-09C9-7347-B37E-CE5CCE65B80A.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0738A4E6-52AC-3142-9CB7-4E059E55D4B4.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/074128AA-4BEC-6C45-9652-F6C93CF64D8E.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/08686EAC-2BF0-9848-8161-F926A417C735.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/087313BD-CAC2-DA4B-B0E8-65B3EC9C100D.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0886CEBB-2006-3140-AFE1-722ACD2E79F1.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/08B336D5-3BD6-CB43-9B87-311488546C76.root',

        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/09181BBE-2688-BF45-9210-E5272F847646.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/094109BF-B509-2E41-8962-A828AD2956B1.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/09413E06-E33D-4941-BB6A-B38E2D3710B0.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/098AF77C-D5BB-4348-AAFB-E186E1730D28.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/09C67C38-22ED-BB44-A84C-DA1E8AC42EE2.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0A0F0808-5DA0-474F-928A-30D916D9AB1B.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0A529A9C-C6DC-B04E-BC81-A64144461C4F.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0A75B49C-CA2D-B249-8226-0EBE36E6B7F9.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0BA377C2-8CA6-254D-889B-789E6AF31377.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0C9F9B98-3AE4-4F46-A1B6-59BAAC2A0A34.root',

        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0CB669FF-60C0-F84F-8165-6D75030DE776.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0CC3059D-4737-B349-8C93-F8330263F87E.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0D0AF398-45F9-4C47-B639-DE12528362BC.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0D2A5724-E01E-CC4B-A76F-1152B4AE56D9.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0DB2F6E9-9E2F-0B48-BC9A-75620901E60D.root',
        '/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0E7EC6AF-3277-CE45-A609-B7186BC7FB7C.root',

        ),##<<>>fileNames = cms.untracked.vstring
)##<<>>process.source = cms.Source("PoolSource"


## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))#ST
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))#T
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))#LT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500))#US
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(12500))#VS
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(25000))#S
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))#SM
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250000))#M
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000000))#ML
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500000))#L
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))#F

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = options.globalTag  

## Create output file
## Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFileName))
				   
# Make the tree 
process.tree = cms.EDAnalyzer("LLPgammaAnalyzer",
   ## flags
   #hasGenInfo = cms.bool("False"),
   ## additional collections
   ## tracks
   tracks = cms.InputTag("unpackedTracksAndVertices"),
   ## vertices
   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
   ## pfcandidates
   pfcandidates = cms.InputTag("packedPFCandidates"),
   ## rho
   rho = cms.InputTag("fixedGridRhoFastjetAll"), #fixedGridRhoAll
   ## METs
   mets = cms.InputTag("slimmedMETs"),
   ## jets
   #jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
   jets = cms.InputTag("slimmedJets"),
   ## electrons
   electrons = cms.InputTag("slimmedElectrons"),
   ## muons
   muons = cms.InputTag("slimmedMuons"),
   ## photons
   gedPhotons = cms.InputTag("slimmedPhotons"),
   ootPhotons = cms.InputTag("slimmedOOTPhotons"),
   ## ecal recHits
   recHitsEB = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
   recHitsEE = cms.InputTag("reducedEgamma", "reducedEERecHits"),
   ## superclusters
   superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters"),
   ootSuperClusters = cms.InputTag("reducedEgamma", "reducedOOTSuperClusters"),
   ## caloclusters
   caloClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters"),
)##<<>>process.tree = cms.EDAnalyzer("LLPgammaAnalyzer"


# Set up the path
process.tree_step = cms.EndPath(process.tree)
process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.content_step = cms.Path(process.content)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.schedule = cms.Schedule(process.tree_step)
process.options = cms.untracked.PSet()

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.
#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
#process = miniAOD_customizeAllData(process)
# End of customisation functions

