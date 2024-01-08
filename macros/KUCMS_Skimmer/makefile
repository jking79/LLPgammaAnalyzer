##################################################3
#
#	Makefile for Root & RJR
#
#################################################3

RTCFLAGS = $(shell root-config --cflags)
RTGLIBS = $(shell root-config --glibs)

RFCFLAGS    = $(shell restframes-config --cxxflags)
RFGLIBS     = $(shell restframes-config --libs)

#BayesianClustering stuff
#paths relative to installation and dir structure of skimmer code
BCCFLAGS = -I$(CMSSW_BASE)/src/BayesianClustering/include
BCCFLAGS  += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0/include/
BCCFLAGS  += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2/include/
BCCFLAGS  += -I/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/json/3.7.3/include/
BCCFLAGS += -frounding-math
BCGLIBS  = -L$(CMSSW_BASE)/src/BayesianClustering/lib -lBayesianClustering
BCGLIBS   += -Wl,-rpath-link,/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/boost/1.63.0/lib/ -lboost_thread
BCGLIBS  += -L/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2/lib/ -lCGAL -lCGAL_Core


CC=g++

p1=runKUCMSAodSkimmer
p2=phoSelSkimHistMaker
p3=skimMultiHistMaker
p4=crStyleSkimHistMaker
p5=skimSVHistMaker

all: $(p1) 

$(p1): $(p1).cpp
	$(CC) -o $(p1).obj $(p1).cpp $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS) $(BCCFLAGS) $(BCGLIBS)

$(p2): $(p2).cpp
	$(CC) -o $(p2).obj $(p2).cpp $(RTCFLAGS) $(RTGLIBS)

$(p3): $(p3).cpp
	$(CC) -o $(p3).obj $(p3).cpp $(RTCFLAGS) $(RTGLIBS)

$(p4): $(p4).cpp
	$(CC) -o $(p4).obj $(p4).cpp $(RTCFLAGS) $(RTGLIBS)

$(p5): $(p5).cpp
	$(CC) -o $(p5).obj $(p5).cpp $(RTCFLAGS) $(RTGLIBS)
