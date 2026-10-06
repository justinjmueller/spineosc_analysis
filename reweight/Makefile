# Build for the event reweighting executable.
#
#   make            builds ./run_reweight
#   make clean      removes it
#
# ROOT supplies both the compile and link flags through root-config, which must
# be on PATH (the detsys calculator is built the same way, by hand).

CXX      ?= g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra $(shell root-config --cflags)
LDLIBS   := $(shell root-config --libs)

# configuration.h is a prerequisite as well as main.cc: it holds the sample
# table and the reweighting list, so editing it alone must still rebuild.
run_reweight: main.cc configuration.h
	$(CXX) $(CXXFLAGS) main.cc -o $@ $(LDLIBS)

clean:
	rm -f run_reweight

.PHONY: clean
