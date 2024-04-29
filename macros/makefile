CC=g++
#CFLAGS = -pthread -std=c++1z -m64 -pthread -lm -ldl -rdynamic
CFLAGS = -pthread -std=c++17 -m64
#LIBS = -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc
LIBS = -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic
#RTLOC = /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_5/external/slc7_amd64_gcc700/bin/../../../../../../../slc7_amd64_gcc700/lcg/root/
#RTLOC = /cvmfs/cms.cern.ch/slc7_amd64_gcc900/cms/cmssw/CMSSW_11_3_0_pre6/external/slc7_amd64_gcc900/bin/../../../../../../../slc7_amd64_gcc900/lcg/root/
RTLOC = /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_20/external/slc7_amd64_gcc700/bin/../../../../../../../slc7_amd64_gcc700/lcg/root/
#RTVRS = 6.12.07-ogkkac4
#RTVRS = 6.22.08-ljfedo
RTVRS = 6.14.09-pafccj6
INCLUDES = -I$(RTLOC)$(RTVRS)/include
LFLAGS = -L$(RTLOC)$(RTVRS)/lib

RTCFLAGS = $(shell root-config --cflags)
RTGLIBS = $(shell root-config --glibs)

p1=hist_maker_aodsim_test
p2=hist_maker_aodsim
p3=wc_ku_InterCali_aveRecHit_isoot_kucc_fliter_red
p4=wc_ku_plot2dResolution_wGausErr_mf
p5=wc_ku_timefitter_wErr
p6=hist_maker_egammares_diag

all: $(p2) 
#$(p2) $(p3) $(p6) $(p22)
# $(p4) $(p5) $(p7)

$(p1): $(p1).cpp
	$(CC) -o $(p1).obj $(p1).cpp $(RTCFLAGS) $(RTGLIBS)

$(p2): $(p2).cpp
	$(CC) -o $(p2).obj $(p2).cpp $(RTCFLAGS) $(RTGLIBS)

$(p3): $(p3).cpp
	$(CC) -o $(p3).obj $(p3).cpp $(RTCFLAGS) $(RTGLIBS)

$(p4): $(p4).cpp
	$(CC) -o $(p4).obj $(p4).cpp $(RTCFLAGS) $(RTGLIBS)

$(p5): $(p5).cpp
	$(CC) -o $(p5).obj $(p5).cpp $(RTCFLAGS) $(RTGLIBS)

$(p6): $(p6).cpp
	$(CC) -o $(p6).obj $(p6).cpp $(RTCFLAGS) $(RTGLIBS)
