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

filelist = [
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022B-PromptReco-v1_352400-358400_dispho/230105_215226/0000/output_23.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022B-PromptReco-v1_352400-358400_dispho/230105_215226/0000/output_9.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_103.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_37.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_39.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_44.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_450.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_453.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_564.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_737.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_886.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_936.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0000/output_996.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1187.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1268.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1281.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1284.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1296.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1319.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1321.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1325.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1328.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1358.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1430.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1432.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1443.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1446.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1502.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1507.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1508.root',
	'store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1513.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1533.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1541.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1551.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1557.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1559.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1576.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1587.root',
	'store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1588.root',
	'store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1590.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1591.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1595.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1596.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1598.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1604.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1606.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1610.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1611.root',
	'store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1618.root',
	'store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1623.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1624.root',
	'store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1626.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1629.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1631.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1632.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1634.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1636.root',
	'store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1678.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022C-PromptReco-v1_352400-358400_dispho/230105_215236/0001/output_1994.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022D-PromptReco-v1_352400-358400_dispho/230105_215246/0000/output_53.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022D-PromptReco-v1_352400-358400_dispho/230105_215246/0000/output_59.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022D-PromptReco-v1_352400-358400_dispho/230105_215246/0000/output_62.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022D-PromptReco-v1_352400-358400_dispho/230105_215246/0000/output_64.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022D-PromptReco-v1_352400-358400_dispho/230105_215246/0000/output_69.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022D-PromptReco-v2_352400-358400_dispho/230105_215255/0000/output_214.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022D-PromptReco-v2_352400-358400_dispho/230105_215255/0000/output_324.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022D-PromptReco-v2_352400-358400_dispho/230105_215255/0000/output_351.root ',
	'store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022G-PromptReco-v1_362300-362800_dispho/230109_183717/0000/output_452.root',
	'store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022G-PromptReco-v1_362300-362800_dispho/230109_183717/0000/output_457.root',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022G-PromptReco-v1_362300-362800_dispho/230109_183717/0000/output_461.root ',
	'/store/user/jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/gammares_tt_kucc_126_v2b_EGamma_MINIAOD_Run2022G-PromptReco-v1_362300-362800_dispho/230109_183717/0000/output_633.root'
]

mspc = '/store/user/jaking/' 
mgrp = '/store/user/lpcsusylep/jaking/' 
eosll = 'eos root://cmseos.fnal.gov rm ' 

for name in filelist :

	command = eosll+name
	#print( command )
	bashout( command )

