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
#outfilename = 'llpgana_t90L_eigen96t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 85 with slope cut at chiprob 0.8, slope not dirct corrted 
#outfilename = 'llpgana_t91L_eigen96t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 85 with slope cut at chiprob 0.95, slope is dirct corrted
#outfilename = 'llpgana_t92ML_eigen96t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 91, w/slope leta*2.2 based & sphericty 0.7- 0.95 cut
#outfilename = 'llpgana_t93ML_eigen96t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 92 w/ exclude leta < 2 && in slope : ltime*10
#outfilename = 'llpgana_t94ML_eigen96t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 93 removed plot cuts in slope calc w/ ltime*1000
#outfilename = 'llpgana_t95L_eigen96t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 94 adjusted axis + added eta v rotated slope plot
#outfilename = 'llpgana_t96L_eigen96t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 95 add phi and eta range plots 
#outfilename = 'llpgana_t97LM_eigen96t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 96 + 
		# photon,oot,electron seperate + plots ( jet matfched and each collection) && rh collection plots && gen dr match v e ratio plot
outfilename = 'llpgana_t98LM_eigen96t60_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root' # as 97 with seed time for p/o/e collections

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

        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/00B87525-94D1-C741-9B03-00528106D15A.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/00F0369F-B9FA-5643-B0F4-F286B80B30B5.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0173DB9C-3E46-4741-B80B-9671E9E917CA.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/017C1808-933C-0D48-95B5-A8B9A0B16601.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/01C71814-B84D-1D4A-8696-61E404084059.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/020C4D59-9FF4-BA4C-B24C-6B77F9F0CC20.root',
        #'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/027298C8-F7C2-7848-B241-DD6F07F41C2F.root',
        ##'/store/data/Run2018A/JetHT/MINIAOD/UL2018_MiniAODv2-v1/260000/0388A30D-E3A0-2947-8F6A-52B8B3F4C381.root',
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

		# AOD

		#'file:465434B8-0A67-5444-8D8D-7E71F561B585.root',#/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/
		
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/0012E808-8EE9-A34F-8745-6BCC6A7EDDA4.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/04F4CD8F-BD5D-8543-A1E5-DC83CE24A0B4.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/07BB6AC8-A8B8-5843-88F0-4415E597BE3E.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/0D394C8D-554D-9147-A625-5A25499F64B4.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/0E00EB5C-F80A-8E42-97F1-5CFDD8577C2C.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/2065699C-D3F4-AB42-8E28-1B76DECE22CC.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/23BA60FD-DCFA-9148-8EE7-E66E06876205.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/2D55C1B7-0BE6-9C49-818C-4AB2C762BA60.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/311CE5A6-3977-724A-8B13-BCCBF52014AB.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/3ABF9B7B-3A90-504E-A9CF-729A5B893100.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/465434B8-0A67-5444-8D8D-7E71F561B585.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/4D33F08D-B498-934C-BD6B-DE16A396F18B.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/63153178-6A55-6D4C-9124-732642B61948.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/64014718-25AA-8F44-B658-B6CE461492FF.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/657E6CEB-4B64-FA44-BFE4-EE9EB6B494AC.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/7414193C-D095-CB43-9F08-496F005C035B.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/7BA17E0D-90DE-A443-90A9-CAF5C4FA6870.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/8413892A-18DF-F84B-9F7B-466F65743588.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/89E366BD-5DB5-E742-B176-3FD47F3C018D.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/9947E79E-2BD1-2F4F-95DC-ED53908ABEE0.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/9B351595-F86C-3E4C-A2D3-2402F25A1FDA.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/A1E97803-960B-9D48-B995-5C5402D68893.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/A299CC74-5753-7241-8264-50914BB49A1B.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/A3171C6D-3098-CC42-BBF6-46C9010C95D9.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/A5936E35-0B2A-5B47-86CA-168060D36FFB.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/B974D226-40AC-2D47-A690-2B83007D0B59.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/C520D3D9-5565-8F45-B79C-8458F4779333.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/C84AE450-446E-F846-805E-04149ED7F006.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/D15AA535-1FD2-0746-90C5-081D951A868F.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/D742AFF4-F7B3-C941-A2FD-753A42D46582.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/D9198D3F-BD8A-C14C-B186-62F43A6CB82B.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/DD64E510-520F-684B-A6B3-C00597BCCDA5.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/DE301D58-F24B-5C49-801C-30CAA33F9BD4.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/E3C58457-CAB5-0142-8D4A-A65C67BDB4BF.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/E8848903-3A7F-984B-9DCF-36827CC691BC.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/EFD7A614-9C2E-1940-B475-2A071218E3AA.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/F0DC76F2-6EA3-B743-9619-8A1AB6CA5479.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2520000/F7F12CFC-4918-AC42-A080-95C2A515A58F.root',

        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2820000/2031AF8F-3D20-754F-989C-935EF1DB438B.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2820000/5CC6B992-94FB-344E-8F9D-B6ECD89B5E4E.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/2820000/78D3AFDF-40DC-7B4D-A514-676968E4DB20.root',

        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/03B27B2E-6187-314E-BB26-E52E815ED0E7.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/069C7211-7C8C-7B4B-8AB3-1689F90722C4.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/0732B1E5-CFB8-6B40-BCED-614F2CC6A744.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/0FD8F44C-85AC-1443-8B1F-0C4E6B917DCE.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/1656B135-9C19-7646-B40A-9D9AC6857BC4.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/189BFC6F-BE84-C848-815E-57ECECCC632D.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/19430A96-3AE7-8B49-897A-346824D70CAC.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/1B166C2A-57C4-BC41-9451-6066E551CDB7.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/1B79B115-C7F8-E441-A7B1-EF12316079FA.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/1DF934C7-F9F6-7D4D-A8C5-780F303802DF.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/1EC5465F-1F82-564E-85E2-73751643F2A9.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/247456F6-E229-6041-BFDC-AFAD8081C84E.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/26650AC6-C1C3-DA44-9959-98EBF3D6021B.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/29C2527F-9484-2D4C-9A94-AE023C4B5AB6.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/2A459101-653A-F04D-9065-CD49FC3B0890.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/327960A5-6266-D549-8961-509374FA788B.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/33D96F68-4FDD-164E-8965-FA118E6D0875.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/35E14C04-02E7-2C49-BE06-10B4EEF9E572.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/36F20E45-6294-7644-85EF-12F2A2F9E7DF.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/3831C3EC-3BE7-5E49-BF23-3C2F3B080F8B.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/387FEAFF-43D2-BC47-81B8-181D9F033222.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/39D3A2E1-2076-6E4F-8B8E-0EF74104AE66.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/3D19EF03-DD9A-4144-833D-B50BD7D7B7E7.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/3D71DE72-7F78-B04F-B7D4-97A60211D38E.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/3E50CA4E-F829-EB40-B572-56BF054860D4.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/3EF9B3BC-6A44-DF41-9A4B-B69D363FEB3F.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/4080C709-C177-9941-925C-9C8F0CE15990.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/421BF6D3-2C1C-7345-ADE6-7B055F83BBAD.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/4295DD83-145C-7745-A375-B62EF1C74444.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/45BA9525-1614-6A41-B7CE-73B191FB5E72.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/470918CC-1BBA-5840-8C76-3862EBE6B4D8.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/4D63BA75-4A21-464C-AFBD-DFD748D708E6.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/4D6E8F7E-00EC-E64D-AF89-1C8956B7047C.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/52CC470A-D803-D547-9B48-30569A57D30E.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/568D1913-D9B2-5149-87BD-F0BCEC308F68.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/588E3061-104F-DF42-9C35-944210E723A5.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/6091B01D-7EB9-0945-9E3B-3CFD68BEE62B.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/61B8DD6F-D71E-6843-B08D-3872F6FBA3FC.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/6289E405-881E-A840-97B7-85106256C631.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/62E39DA9-4F1C-6647-B123-C2E3BEB7B533.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/6450FC1A-329A-AE47-B017-E34FB733A9AE.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/6758CDE6-BEA5-684A-BDC6-7F963C77D9C8.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/67C8F6A1-1B2B-B74D-AE69-D23052935884.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/696E278D-29AD-7E46-B6E0-F26AB8A075D9.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/6BF30361-6A8E-0347-A72C-9040F4B4E262.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/6D250534-3894-6D48-A22F-80DAFA5D56CD.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/6DB45F6F-DB51-474D-AC9D-2FC4105DD9BB.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/6E478999-2F63-3B40-B9F5-56FF650C05BF.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/70C13508-56B6-5B4E-AB23-717A3F1084DA.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/79E44F9C-EF42-7D4F-BBC1-860C2075AB92.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/7B93749A-427C-C34A-8152-1DBE8A6B58AC.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/7EB8F765-D2FE-7C4D-8DDE-34EC1C2C8409.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/96A0074F-A767-7B46-925B-AD40AD82BA4C.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/96F795EF-6EB3-4C4A-90E8-8653330ADC4D.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/98445316-AA4A-F942-B503-A939B0DAF9B8.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/98771765-BDC7-B24D-91F2-367A220B0FE2.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/A064F863-5647-7F4F-8A0D-B8F65D1B1120.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/A15E834F-3723-2641-B5B2-D4F29F6CA21A.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/A4C87AAA-CF23-7C47-8042-D8133A13B39A.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/A6CF8B46-F089-C647-9278-B203C22E90A6.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/AD4B8B4A-B89B-CB4A-9C4F-90ECA2E55D64.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/B0E466AF-2688-2041-918F-9D592575B219.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/B770BFCC-1049-A446-808B-E9DC42AC9A6D.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/B783DFCD-92E9-C343-B760-B482383DECED.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/B836DF57-CEA1-0B4E-B29C-CD0B7CCA4199.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/B87FE964-3217-4645-89F9-DEB68879DCB6.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/B909156B-55D7-2B4F-8150-631AB8197681.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/BB24C660-C5BC-D840-A581-A345D7EB7A36.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/BCBCAF43-810D-4748-8870-7B746B1F3CF1.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/BCDABA9C-D9FB-7047-A1CF-5403C46803ED.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/CFBAE73E-DADA-1D47-A8EC-C7E22243FA1D.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/D5A00180-E4EB-9741-9BEA-64474DB1E08D.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/D926000B-B22E-5F46-8E6D-81839745E4FF.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/DC18C575-128E-F143-A79A-63542C02BB86.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/DEE842F8-4477-2848-A37C-ECA0FF5FA43E.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/E275542E-B954-DF43-B7AF-6144ACF8B60D.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/EAFA49C7-F972-9B43-B39F-606F6B9B79ED.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/EE87B456-BA4E-4D45-8DEA-EBCB73BFEFCE.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/F035FDA0-1275-2C45-994A-E23706179676.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/F3D1B877-A897-444F-947D-FEBC2CD1D710.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/F5144E54-6557-3B4C-90A4-B94FD4D7D930.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/F5583974-DAAC-6141-8E35-805CF88E82D6.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/F65EFF8A-9EC4-A545-9E6E-5B4DEACBBAAB.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/F8992F1A-B01E-A341-BCE2-B67EDA9735A8.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/FB297DBF-77CD-D64F-BF9D-2693C8EB9D59.root',
        '/store/data/Run2018B/JetHT/AOD/15Feb2022_UL2018-v1/28210000/FD318109-957F-9246-8942-877D2B44F106.root',

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
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250000))#M
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000000))#LM
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1750000))#L
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500000))#VL
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

