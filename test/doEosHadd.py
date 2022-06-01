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

mspc = '/store/user/jaking/'
eosll = 'eos root://cmseos.fnal.gov ls '
command = eosll+mspc+'LLPGamma/'
version = '_v13_'
#command = 'ls'

targdirs = []
subdirlist1 = []
subdirlist2 = []
theFileList = ''
dirls = bashout( command ).splitlines()
#print( dirls )
#print( '-------------------------------------------------')
dirselect = 'HTo2LongLivedTo4b'
for line in dirls:
	#print( line )
	if dirselect in line :
		targdirs.append( line )
#print( targdirs )
for mydir in targdirs:
	command1 = command+mydir+'/'
	subdir1 = bashout( command1 ).rstrip().splitlines()
	#print( subdir1 )
	#print( mydir+'/'+subdir1+'/' )
	for line in subdir1 : 
		#print( line )
		if version in line : subdirlist1.append( mydir+'/'+line+'/' )
	#print( subdirlist1 )
	for thesubdir in subdirlist1 :
		command2 = command+thesubdir
		subdir2 = bashout( command2 ).rstrip()
		#print( thesubdir+subdir2+'/0000/' )
		subdirlist2.append(thesubdir+subdir2+'/0000/')
		#for subdir2 in subdirlist2:
		#	command3 = command+subdir2
		#	files = bashout( command3 ).splitlines()
		#	for thefile in files : 
		#		#print( subdir2+thefile )
		#		theFileList +=(subdir2+thefile)

#print( theFileList )
for subdir2 in subdirlist2:
	filename = 'tmp_'+subdir2.split('/')[1]+'.root '
	#print( filename )
	#lists = bashout( "eosls "+mspc+"LLPGamma/"+subdir2 ).rstrip()
	#print( subdir2 )
	haddcommand = "hadd -f "+filename+"`xrdfs root://cmseos.fnal.gov ls -u "+mspc+"LLPGamma/"+subdir2+" | grep '\.root'`"
	#haddcommand = "`xrdfs root://cmseos.fnal.gov ls -u "+mspc+"LLPGamma/"+subdir2+" | grep '\.root'`"
	#haddcommand = "`xrdfs root://cmseos.fnal.gov ls -u "+mspc+"LLPGamma/"+subdir2+"`"
	#haddcommand = "hadd "+filename+lists
	print( haddcommand )
	doCommand( haddcommand )
	print( '---------------------------------------------------' )



	
#print( bashout( 'hadd llpgana_HTo2LongLivedTo4b_t37MC_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root ' + theFileList ) )	
