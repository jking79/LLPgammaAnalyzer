#! /usr/bin/env python

import os
from optparse import OptionParser

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = 'submit',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = 'myworkingArea',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def docrab( dataset ):

    options = getOptions()

    # The submit command needs special treatment.
    if options.crabCmd == 'submit':

        # External files needed by CRAB
        #inputJSON    = 'golden2016.json'
        #inputJSON    = 'golden2017.json'
        inputJSON    = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        from CRABClient.UserUtilities import config
        config = config()

        config.General.workArea    = options.workArea
        config.General.requestName = None

        config.JobType.pluginName  = 'Analysis'
        config.JobType.psetName    = 'llpgana_mc.py'
        config.JobType.pyCfgParams = None

        config.Data.inputDataset   = None
        #config.Data.lumiMask       = inputJSON    # Comment out for MC only
        #config.Data.splitting     = 'Automatic'
        #config.Data.unitsPerJob   = 2000
        config.Data.splitting    = 'EventAwareLumiBased' # MC
        config.Data.unitsPerJob  =  10000 # MC

        config.JobType.allowUndistributedCMSSW = True
        config.Data.publication    = False
        #config.Site.storageSite    = 'T2_US_Nebraska'
        config.Site.storageSite    = 'T3_US_FNALLPC'
        config.Data.outLFNDirBase  = '/store/user/jaking/LLPGamma/'
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDataAndOpts = [[dataset[0]]]

        for inDO in inputDataAndOpts:
            print( '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>' )
            print( 'Input dataset for Crab Job : ' )
            #print( inDO )
            # inDO[0] is of the form /A/B/C. Since A+B is unique for each inDS, use this in the CRAB request name.
            primaryDataset = (inDO[0].split('/')[1]).split('_T')[0]
            print( primaryDataset )
            runEra         = ((inDO[0].split('/')[2]).split('_')[0]+'_'+(inDO[0].split('/')[2]).split('_')[1]).split('MiniAOD')[0]
            print( runEra )
            dataset        = inDO[0].split('/')[3]
            print( dataset )

            trial          = "llpga_v1" # 4 Feb 22 : t37L_ phsc & elesc _005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2 

            config.General.requestName   = trial+"_"+primaryDataset+"_"+runEra+"_"+dataset+"_llpga"
            config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_llpga"

#>>>>>>>>>>>>>>>>>>>     #2018   #globalTag=106X_dataRun2_v28
            #config.JobType.pyCfgParams   = ['globalTag=106X_dataRun2_v28','outputFileName=output.root']
#>>>>>>>>>>>>>>>>>>>	    #2017   #globalTag=106X_dataRun2_v20
            #config.JobType.pyCfgParams   = ['globalTag=106X_dataRun2_v20', 'outputFileName=output.root']
#>>>>>>>>>>>>>>>>>>>	    #2016  #globalTag=106X_dataRun2_v27
            #config.JobType.pyCfgParams   = ['globalTag=106X_dataRun2_v27', 'outputFileName=output.root']
#>>>>>>>>>>>>>>>>>>>     #MC   #globalTag=112X_mcRun3_2021_realistic_v16  #  <<<<<<<   comment/uncomment lumi mask when using/!using MC  >>>>>>>>>>>>>
            config.JobType.pyCfgParams   = ['globalTag=112X_mcRun3_2021_realistic_v16','outputFileName=output.root']

            config.Data.inputDataset     = inDO[0]
            # Submit.
            try:
                print "Submitting for input dataset %s" % primaryDataset + '_' + runEra + '_' + dataset
                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
                os.system("rm -rf %s/crab_%s/inputs" % (config.General.workArea, config.General.requestName))
            except HTTPException as hte:
                print "Submission for input dataset %s failed: %s" % (inDO[0], hte.headers)
            except ClientException as cle:
                print "Submission for input dataset %s failed: %s" % (inDO[0], cle)

    # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print "-"*len(msg)
            print msg
            print "-"*len(msg)
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers)
            except ClientException as cle:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle)


##33333333333333333333333333333333333333333333333333333333333

def run_multi():
 
        run = 'RunIISummer19UL18RECO-106X_upgrade2018_realistic_v11' 
        tune = '_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8' 
        mcdset = 'AODSIM' 

        datasets = [

# Dataset: JetHT UL2018

			#['/JetHT/Run2018A-UL2018_MiniAODv2-v1/MINIAOD',''],
			#['/JetHT/Run2018B-UL2018_MiniAODv2-v1/MINIAOD',''],
			#['/JetHT/Run2018C-UL2018_MiniAODv2-v1/MINIAOD',''],
			#['/JetHT/Run2018D-UL2018_MiniAODv2-v1/MINIAOD',''],
			#['/JetHT/Run2018D-UL2018_MiniAODv2-v2/MINIAOD',''],


# Dataset: QCD HT MC 

            #['/QCD_HT50to100'   +tune+'/'+run+'_L1v1-v2/AODSIM',''], 
            #['/QCD_HT100to200'  +tune+'/'+run+'_L1v1-v2/AODSIM',''],
            #['/QCD_HT200to300'  +tune+'/'+run+'_L1v1-v2/AODSIM',''], 
            #['/QCD_HT200to300'  +tune+'/'+run+'_L1v1_ext1-v1/AODSIM',''], 
            #['/QCD_HT300to500'  +tune+'/'+run+'_L1v1-v2/AODSIM',''], 
            #['/QCD_HT500to700'  +tune+'/'+run+'_L1v1-v2/AODSIM',''], 
            #['/QCD_HT700to1000' +tune+'/'+run+'_L1v1-v2/AODSIM',''], 
            #['/QCD_HT1000to1500'+tune+'/'+run+'_L1v1-v2/AODSIM',''], 
            #['/QCD_HT1500to2000'+tune+'/'+run+'_L1v1-v2/AODSIM',''],
            #['/QCD_HT2000toInf' +tune+'/'+run+'_L1v1-v2/AODSIM',''],


# Dataset: /EGamma/Run2018-12Nov2019_UL2018-/MINIAOD

            #['/EGamma/Run2018A-12Nov2019_UL2018-v2/MINIAOD',''],
            #['/EGamma/Run2018B-12Nov2019_UL2018-v2/MINIAOD',''],
            #['/EGamma/Run2018C-12Nov2019_UL2018-v2/MINIAOD',''],
            #['/EGamma/Run2018D-12Nov2019_UL2018-v4/MINIAOD',''],
    
# Dataset: /DoubleEG/Run2016-21Feb2020_UL2016-/MINIAOD

            #['/DoubleEG/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD',''],
            #['/DoubleEG/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD',''],
            #['/DoubleEG/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD',''],
            #['/DoubleEG/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD',''],
            #['/DoubleEG/Run2016F-21Feb2020_UL2016-v1/MINIAOD',''],
            #['/DoubleEG/Run2016G-21Feb2020_UL2016-v1/MINIAOD',''],
            #['/DoubleEG/Run2016H-21Feb2020_UL2016-v1/MINIAOD',''],

# Dataset: /DoubleEG/Run2017-09Aug2019_UL2017-/MINIAOD

            #['/DoubleEG/Run2017B-09Aug2019_UL2017-v1/MINIAOD',''],
            #['/DoubleEG/Run2017C-09Aug2019_UL2017-v1/MINIAOD',''],
            #['/DoubleEG/Run2017D-09Aug2019_UL2017-v1/MINIAOD',''],
            #['/DoubleEG/Run2017E-09Aug2019_UL2017-v1/MINIAOD',''],
            #['/DoubleEG/Run2017F-09Aug2019_UL2017-v1/MINIAOD',''],

# Dataset: /SingleElectron/Run2017-09Aug2019_UL2017-/MINIAOD

             #['/SingleElectron/Run2017B-09Aug2019_UL2017-v1/MINIAOD',''],
             #['/SingleElectron/Run2017C-09Aug2019_UL2017-v1/MINIAOD',''],
             #['/SingleElectron/Run2017C-09Aug2019_UL2017_EcalRecovery-v1/MINIAOD',''],
             #['/SingleElectron/Run2017D-09Aug2019_UL2017-v1/MINIAOD',''],
             #['/SingleElectron/Run2017E-09Aug2019_UL2017-v1/MINIAOD',''],
             #['/SingleElectron/Run2017F-09Aug2019_UL2017_rsb-v2/MINIAOD',''],
             #['/SingleElectron/Run2017F-09Aug2019_UL2017_EcalRecovery-v1/MINIAOD',''],

            ]

        datasetsMC = [

# Dataset: HTo2LongLivedTo4b

            ['/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-9000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-15000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-1500mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-30000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-10000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-10000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-1000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-250_MFF-60_CTau-500mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-10000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-1000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-10000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-1000mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'],
            ['/HTo2LongLivedTo4b_MH-350_MFF-80_CTau-500mm_TuneCP5_14TeV-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM']

            ]

#	    for dataset in datasets :
        for dataset in datasetsMC :
		    docrab( dataset )



run_multi()

