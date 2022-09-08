## Makefile for community detection code in bipartite networks

CC = gcc
CXX = g++

OPTFLAGS = -Ofast
CFLAGS = $(OPTFLAGS)
CXXFLAGS = $(CFLAGS) -std=c++11 -DUSE_32_BIT_GRAPH
LDFLAGS="-L/opt/homebrew/opt/llvm/lib -Wl,-rpath,/opt/homebrew/opt/llvm/lib"

GOBJFILES = Main.o Timer.o Graph.o Node.o MetaNode.o Community.o biLouvainMethod.o biLouvainMethodMurataPN.o FuseMethod.o


GTARGET = biLouvain

all: $(GTARGET)

$(GTARGET):  $(GOBJFILES)
	$(CXX) $^ $(OPTFLAGS) -o $@

$(FTARGET):  $(FOBJFILES)
	$(CXX) $^ $(OPTFLAGS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

.PHONY: clean

clean:
	rm -f *~ $(GOBJFILES) $(GTARGET) $(FTARGET) $(FOBKFILES)
