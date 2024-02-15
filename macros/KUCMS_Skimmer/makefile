##################################################3
#
#	Makefile for Root & RJR
#
#################################################3

RTCFLAGS = $(shell root-config --cflags)
RTGLIBS = $(shell root-config --glibs)

RFCFLAGS    = $(shell restframes-config --cxxflags)
RFGLIBS     = $(shell restframes-config --libs)

CC=g++

p1=runKUCMSAodSkimmer
p2=phoSelSkimHistMaker
p3=skimMultiHistMaker
p4=crStyleSkimHistMaker
p5=skimSVHistMaker
p6=calibration_averhtime_isoot_filter
p7=stackSkimHistMaker

all: $(p1) 

$(p1): $(p1).cpp
	$(CC) -o $(p1).obj $(p1).cpp $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS)

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

$(p7): $(p7).cpp
	$(CC) -o $(p7).obj $(p7).cpp $(RTCFLAGS) $(RTGLIBS)
