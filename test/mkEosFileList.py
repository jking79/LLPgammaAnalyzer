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

debug = True
#debug = False

mspc = '/store/user/jaking/'
mdis = '/store/user/lpcsusylep/jaking/ecalTiming/EGamma/'
edis = '/eos/uscms/store/mc/RunIIFall17DRPremix/'
eosll = 'eos root://cmseos.fnal.gov ls '
#command = eosll+mspc+'LLPGamma/'
#command = eosll+mspc+'ecalTiming/'
command = eosll+edis
#version = '_v11_'
#version = '_noOOTAmp_'
#version = '_wthOOTAmp_'
version = ''
folder = ''
rootfile = '.root'
dirselect = 'GMSB'
#dirselect = 'HTo2LongLivedTo4b'
#dirselect = '_jitdif_'

targdirs = []
subdirlist1 = []
subdirlist2 = []
subdirlist3 = []
filelist = []
theFileList = ''

dirls = bashout( command ).splitlines()
#print( '-------------------------------------------------')
for line in dirls:
	if dirselect in line : targdirs.append( line )
for mydir in targdirs:
	command1 = command+mydir+'/'
	subdir1 = bashout( command1 ).rstrip().splitlines()
	for line in subdir1 : 
		if version in line : subdirlist1.append( mydir+'/'+line+'/' )
for thesubdir in subdirlist1 :
	command2 = command+thesubdir
	subdir2 = bashout( command2 ).rstrip().splitlines()
	for subdir in subdir2 : subdirlist2.append(thesubdir+subdir+'/')
for subdir2 in subdirlist2:
	lists = bashout( command+subdir2 ).rstrip().splitlines()
	for line in lists :
		if folder in line : subdirlist3.append(subdir2+line+'/')
for subdir3 in subdirlist3:
    lists = bashout( command+subdir3 ).rstrip().splitlines()
    for line in lists :
        if rootfile in line : filelist.append(subdir3+line)
for thefile in filelist:
	print( thefile )

