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
options.register('outputFileName','llpgana_output_t12_id3pt200rh32eta15.root',VarParsing.multiplicity.singleton,VarParsing.varType.string,'output file name created by cmsRun');

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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

## Define the input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#'file:jwk_reco_data_DIGI2RAW.root'),
        #'/store/data/Run2018D/EGamma/MINIAOD/22Jan2019-v2/70001/F0FE59FC-F29E-904B-A1BC-817C1CB09A7E.root'
        '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/062A29F3-416F-6A43-8E0A-90BE80607677.root',
        '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/0632DF20-D859-EE4B-8FC4-3ECC56B2987D.root',
        '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/0685BE81-54CF-3A4C-8D66-8E0E059D263F.root',
        '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/06CBEAC6-F3C1-BE48-96D9-A7EF3CF42048.root',
        '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/06CC5637-77D3-8145-B3BD-0AF5F16410AB.root',
        '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/06DE38F0-4330-C546-8D1A-0865C50FAA14.root',
        '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/071F2728-A0CA-9446-A1E3-38B4FB08D6AE.root',
        '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/077E7C54-C588-F642-BBA8-0960D5FD7A96.root',
        '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/07971C26-63F0-D943-BA77-88D6F636F780.root',
        '/store/data/Run2018A/EGamma/MINIAOD/12Nov2019_UL2018-v2/100000/07A40762-F4E1-9E4D-8750-7DD65EA01B8A.root',
        ),
)


## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = options.globalTag  

## Create output file
## Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", 
		                   fileName = cms.string(options.outputFileName))
				   
# Make the tree 
process.tree = cms.EDAnalyzer("LLPgammaAnalyzer",
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
)


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

