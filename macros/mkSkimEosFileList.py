import subprocess
import sys
import os

# uses python3

def bash( bashCommand ):
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	#process = subprocess.Popen(bashCommand.split())
	output, error = process.communicate()
	return output ,error

def bashout( command ):
	#output = subprocess.check_output( command, shell=True)
	output = subprocess.run( command, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True )
	return output.stdout	

def doCommand( command ):
	output = os.system( command )
	return output

mspc = '/store/user/jaking/'
mdis = '/store/user/lpcsusylep/jaking/'
eosll = 'eos root://cmseos.fnal.gov ls '
#command = eosll+mspc+'LLPGamma/llpga_GMSB_AOD_v48/'
#command = eosll+mspc+'A/'
#command = eosll+mdis+'LLPGamma/llpga_GMSB_AOD_v58/'
#command = eosll+mdis+'/ecalTiming/tt_KUCCRes_126_Test/EGamma/'
#command = eosll+mdis+'LLPGamma/llpga_GJets_AOD_v57/'
#command = eosll+mdis+'LLPGamma/llpga_GMSB_AOD_v59/'
#command = eosll+mspc+'ecalTiming/gammares_ttcc_131_v11_diag/'
#command = eosll+mdis+'ecalTiming/EGamma/'
#command = eosll+mspc+'EGamma/'
#command = eosll+mdis+'KUCMSNtuple/GMSB_AOD_v1/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GJETS_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_WJETS_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_QCD_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GMSB_AOD_v6/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_ZJETS_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_DYTT_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_JetHTi_Met50_AOD_v2/'
#command = eosll+mdis+'EcalTiming/ZeroBias/'
command = eosll+mspc+'ecalTiming/EGamma0/'

version = ''
#version = '_v11_'
#version = '_noOOTAmp_'
#version = '_wthOOTAmp_'
rootfile = '.root'
#dirselect = 'HTo2LongLivedTo4b'
#dirselect = '_newRtParams4_v26b_'
#dirselect = '_newRtParams3_test_v26_'
#dirselect = 'tt_kurhs_124cc5_cert'
#dirselect = '22eraC_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = 'noOOTCC_kustc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = '22eraC_CCstc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = 'noOOTCC_kustc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357101-357268'
#dirselect = 'CCstc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357101-357268'

#dirselect = 'GMSB'
#dirselect = 'AOD'
#dirselect = 'WJetsToLNu_HT-800'
#dirselect = 'QCD_HT100to200'
#dirselect = 'GMSB_L-400TeV'
#dirselect = 'DYJetsToLL_M-50'
#dirselect = 'TTJets'
dirselect = 'v2_EGamma'
#dirselect = 'ecaltiming_dqm_132r3prompt3'

#dirselect = ''

#debug = True
debug = False

#deep = True
deep = False

targdirs = []

dirls = bashout( command ).splitlines()
print( '************************************************')
for line in dirls:
	#print( line )
	if dirselect in line : targdirs.append( line )
    #targdirs.append( line )
print( targdirs )

for line2 in targdirs :

    subdirlist1 = []
    subdirlist2 = []
    subdirlist3 = []
    filelist = []
    theFileList = ''

    #for mydir in targdirs:
    print( '-------------------------------------------------')
    print( line2 )
    print( '-------------------------------------------------')
    #subdirlist1.append( line+'/' )
    
    #if debug : print( subdirlist1 )
    #for thesubdir in subdirlist1 :
    thesubdir = line2
    command2 = command+thesubdir+'/'
    if debug : print( command2 )
    subdir2 = bashout( command2 ).rstrip().splitlines()
    for subdir in subdir2 : 
    	command3 = command+thesubdir+'/'+subdir+'/'
    	subdir3 = bashout( command3 ).rstrip().splitlines()
    	for subsubdir in subdir3 : 
    		subdirlist2.append(thesubdir+'/'+subdir+'/'+subsubdir)
   
    if debug : print( subdirlist2 ) 
    for thesubdir2 in subdirlist2 :
        command4 = command+thesubdir2+'/'
        subdir4 = bashout( command4 ).rstrip().splitlines()
        #print( thesubdir+subdir2+'/0000/' )
        for subdir in subdir4 :
            subdirlist3.append(thesubdir2+'/'+subdir+'/')
            #command5 = command+thesubdir+subdir+'/'
            #subdir5 = bashout( command5 ).rstrip().splitlines()
            #for subsubdir in subdir5 :
                #subdirlist3.append(thesubdir+subdir+'/'+subsubdir+'/')

    if debug : print( subdirlist3 )
    for subdir2 in subdirlist3:
    	lists = bashout( command+subdir2 ).rstrip().splitlines()
    	for lline in lists :
    		if rootfile in lline : filelist.append(subdir2+lline)

    select =  line2.split("Tune")
    outfile = 'kuntuple_' + select[0] + '_v6.txt'
    print( outfile )
    outf = open( outfile, 'w' )
    filelist = subdirlist3
    for thefile in filelist:
    	outf.write( thefile + '\n' )
    outf.close()



