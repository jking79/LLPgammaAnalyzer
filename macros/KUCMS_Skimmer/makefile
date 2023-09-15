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
p2=baseSkimHistMaker

all: $(p1) 

$(p1): $(p1).cpp
	$(CC) -o $(p1).obj $(p1).cpp $(RTCFLAGS) $(RTGLIBS) $(RFCFLAGS) $(RFGLIBS)

$(p2): $(p2).cpp
	$(CC) -o $(p2).obj $(p2).cpp $(RTCFLAGS) $(RTGLIBS)
