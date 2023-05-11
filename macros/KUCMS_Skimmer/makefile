##################################################3
#
#	Makefile for Root & RJR
#
#################################################3

INCLUDES = $(shell root-config --cflags)
LFLAGS = $(shell root-config --glibs)

CC=g++

p1=runKUCMSAodSkimmer


all: $(p1) 

$(p1): $(p1).cpp
	$(CC) -o $(p1).obj $(p1).cpp $(INCLUDES) $(LFLAGS)
