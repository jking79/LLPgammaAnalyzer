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
#options.register('globalTag','106X_dataRun2_v28',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','112X_mcRun3_2021_realistic_v16',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
options.register('globalTag','94X_mc2017_realistic_v14',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#94X_mc2017_realistic_v14
#102X_mc2017_realistic_v7
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
#outfilename = 'llpgana_GMSB_t80M_pheigen95t60r9_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 79 + genJetDrmatch
outfilename = 'llpgana_mc_MINIAODSIM_GMSB100TeV_v80a_2per_FL.root'

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
#process.MessageLogger.destinations = ['cout', 'cerr']
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 10000

## Define the input source

GMSBpath = '/store/mc/RunIIFall17MiniAODv2/' 
gmsbver = 'MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/'
gmsbver2 = 'MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/' 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

		#'file:jwk_reco_data_DIGI2RAW.root'),

         ## HTo2LongLivedTo4b

        #'/store/mc/Run3Winter21DRMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_14TeV-pythia8/MINIAODSIM/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280000/17bd2d71-8a76-46c5-947a-7ea2b1df44b6.root'

		 ## GMSB

        #'/store/mc/RunIIFall17MiniAODv2/GMSB_L-150TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/60000/B03A864D-79BA-E811-BA05-0242AC1C0500.root'

        GMSBpath+'GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/'+gmsbver2+'10000/841F7B7B-F629-E911-88C8-0025905B860C.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/'+gmsbver2+'10000/3C4B6292-F629-E911-8AD4-24BE05CECBD1.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/'+gmsbver2+'100000/D428E60B-AD32-E911-B3EA-0025905C42A8.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/'+gmsbver2+'100000/6EE235DE-5433-E911-B06D-0CC47AFF025C.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/'+gmsbver2+'100000/FA2A3E25-0032-E911-AEFB-00266CF8586C.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/'+gmsbver2+'100000/BEF92BA4-5433-E911-BDEB-44A842CFC9E6.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/'+gmsbver2+'100000/AA8BA6A5-5533-E911-9468-B083FED42FAF.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/'+gmsbver2+'100000/94C0303E-5533-E911-851A-24BE05C6E7E1.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/'+gmsbver2+'100000/BCDCE0E5-5433-E911-AB00-001E67F33631.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/'+gmsbver2+'100000/2E56CD8A-5533-E911-8E63-00259019A41E.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/'+gmsbver+'270000/0C22794B-F7DA-E811-97AF-001EC94B4F57.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/'+gmsbver+'270000/4C422BFD-F6DA-E811-A790-001E67DDCC81.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/'+gmsbver+'270000/BC90DCFB-F6DA-E811-8186-00259075D62E.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/'+gmsbver+'270000/C67BE728-F7DA-E811-9E36-B499BAAC0572.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/026A1D51-6DB4-E811-90E0-0025905D1CB4.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/5A6BB0E2-6CB4-E811-A4B2-E0071B73C630.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/D0579E44-3DB4-E811-A58E-1C6A7A26E465.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/D4AC52BB-6DB4-E811-8F2D-A4BF010120C5.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/6E9F4FA2-67AE-E811-9E39-24BE05CEECD1.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/FA10588B-55B3-E811-A08D-A4BF0101DDD7.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/D6DC264C-55B3-E811-9607-001E67E63AE6.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/10191D17-55B3-E811-B24A-EC0D9A82261E.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/56778B4A-55B3-E811-A0E7-0CC47AF9B32A.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/6E777363-55B3-E811-98A0-0CC47A13D416.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/BEEB6449-55B3-E811-8F9D-A0369FE2C0D2.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/DED58E6F-55B3-E811-982A-008CFA0087C4.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/'+gmsbver+'120000/BA26173B-BCBC-E811-9C4B-FA163EC2B778.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/'+gmsbver+'120000/FA59BD8A-2FC7-E811-A4C7-008CFA197E18.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/C25F11E2-A1B3-E811-B4A4-0090FAA58924.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/EAA4D0DD-A1B3-E811-A44A-00269E95B0A4.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/304BEDCA-A1B3-E811-8A43-0023AEEEB208.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/7A5A6CE8-A1B3-E811-989F-B499BAAC0626.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/BC34C0E1-A1B3-E811-8A82-0CC47AF9B1AA.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/4E1FFED6-A1B3-E811-BC6B-002590725380.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/A6DCC1C0-43B3-E811-B2FE-24BE05C48801.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/'+gmsbver+'120000/C4992F84-16D9-E811-9313-D4AE526A10E8.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/'+gmsbver+'120000/6806C678-16D9-E811-ADC1-AC1F6B8DD22E.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/16678BE8-90B6-E811-BE66-484D7E8DF085.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/7A3167FA-88B6-E811-AECA-AC1F6B0DE454.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/827D4CB4-8BB6-E811-B294-246E96D10CBC.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/E67A2621-89B6-E811-8183-24B6FDFFBC48.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/80ED4776-C2B3-E811-8454-E0071B6CAD00.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/DAF3347D-9DB6-E811-885F-90B11CBCFF4E.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/AED6BBBC-9DB6-E811-81C0-FA163E3C0C41.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/1E320180-2AB0-E811-B097-842B2B76653D.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/F6F9A29D-2AB0-E811-8A2F-A4BF01159514.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/8C02DA9A-2AB0-E811-A4B0-0025905C5502.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/0814AB87-2AB0-E811-82C7-E0071B7A7830.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/E25298F9-2AB0-E811-A54C-009C02AABB60.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/FC0993AB-2AB0-E811-AE73-0025905B8564.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/CA16B7E8-2AB0-E811-93A6-1C6A7A21A6B3.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/'+gmsbver+'60000/26DE9DF3-BBD7-E811-A0C6-A4BF0101DB93.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/'+gmsbver+'60000/A2F6CF9B-BBD7-E811-A94D-0425C5DE7BF2.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/B4FF3ECF-52B5-E811-A6D9-0025905C3DF8.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/9A1B0A02-53B5-E811-A72F-001F29089F2A.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/4E1493A8-53B5-E811-8B2D-5065F3812261.root',
        GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/4CDE3601-53B5-E811-A342-1866DA890B94.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/B498C3D9-52B5-E811-994A-001E67E6F8FA.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/202BC5D4-52B5-E811-827C-0CC47AD98D00.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/A0F1EAFE-52B5-E811-9D20-AC1F6B8DBE36.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'80000/58EAF902-53B5-E811-82C9-0026B9278651.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'270000/C456A226-0DD9-E811-BB19-AC1F6B1AF142.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'270000/3CD4E9BE-F6D8-E811-891A-90E2BAD57CD0.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'270000/3E6187BA-F6D8-E811-ACF5-6CC2173D6B10.root',
        #GMSBpath+'GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/'+gmsbver+'100000/6C2B5185-0BD8-E811-80B8-0CC47A1E0484.root',


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
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))#ST
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))#T
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))#LT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500))#US
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(12500))#VS
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(25000))#S
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))#SM
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250000))#M
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500000))#L
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))#F

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = options.globalTag  

## Create output file
## Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", 
		                   fileName = cms.string(options.outputFileName))
				   
# Make the tree 
process.tree = cms.EDAnalyzer("LLPgammaAnalyzer",
   ## flags
   hasGenInfo = cms.bool(options.hasGenInfo),
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

process.options = cms.untracked.PSet( 
    #numberOfThreads = cms.untracked.uint32(4), 
	#numberOfStreams = cms.untracked.uint32(4), 
    SkipEvent = cms.untracked.vstring('ProductNotFound'), 
    #wantSummary = cms.untracked.bool(True) 
)#process.options = cms.untracked.PSet

#do not add changes to your config after this point (unless you know what you are doing)
#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#process=convertToUnscheduled(process)

# customisation of the process.
#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
#process = miniAOD_customizeAllData(process)
# End of customisation functions

