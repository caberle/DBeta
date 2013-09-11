#Make File for DBeta Generator
CXXFLAGS += -I/home/christoph/clhep-2.0.4.5/installation_dir/include
CXXFLAGS += $(shell root-config --cflags)
ROOTLIB := $(shell root-config --libs)
LIBS += -L. -L/home/christoph/clhep-2.0.4.5/installation_dir/lib -lCLHEP-2.0.4.5 $(ROOTLIB)

%:%.cc
	$(CXX) $(CXXFLAGS) $(ROOTCFLAG) $< -o $@ $(LIBS)

TARGETS= DBeta

all: $(TARGETS)

clean:
	rm -rf $(TARGETS)
