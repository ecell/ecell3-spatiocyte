OBJS=\
SpatiocyteVisualizer

DMS=\
SpatiocyteStepper\
SpatiocyteProcess\
MassActionProcess\
ReactionProcess\
VisualizationLogProcess\
IteratingLogProcess\
DiffusionInfluencedReactionProcess\
OscillationAnalysisProcess\
SpatiocyteNextReactionProcess\
FluorescentImagingProcess\
CoordinateLogProcess\
MoleculePopulateProcess\
DiffusionProcess\
PolymerizationParameterProcess\
PolymerizationProcess\
PeriodicBoundaryDiffusionProcess\
PolymerFragmentationProcess

ECELL3_DMC = ecell3-dmc
CXX = g++
CXXFLAGS = -Wall -O3 -g
CXXFLAGS += $(shell pkg-config --cflags gtkmm-2.4 gtkglextmm-x11-1.2)
CPPFLAGS = -DG_DISABLE_DEPRECATED -DGDK_PIXBUF_DISABLE_DEPRECATED -DPNG_SKIP_SETJMP_CHECK # -DGDK_DISABLE_DEPRECATED 
GUILIBS =
GUILIBS += $(shell pkg-config --libs gtkmm-2.4 gtkglextmm-x11-1.2)
SPATIOCYTE = spatiocyte
OBJECTS=${OBJS:=.o}
SOS=${DMS:=.so}

all:	$(SOS) $(SPATIOCYTE)

VisualizationLogProcess.so: 	VisualizationLogProcess.cpp
	$(ECELL3_DMC) -o VisualizationLogProcess.so --ldflags=SpatiocyteProcess.so VisualizationLogProcess.cpp

FluorescentImagingProcess.so: 	FluorescentImagingProcess.cpp
	$(ECELL3_DMC) -o FluorescentImagingProcess.so --ldflags="SpatiocyteProcess.so VisualizationLogProcess.so" FluorescentImagingProcess.cpp

ReactionProcess.so: 	ReactionProcess.cpp
	$(ECELL3_DMC) -o ReactionProcess.so --ldflags=SpatiocyteProcess.so ReactionProcess.cpp

IteratingLogProcess.so: 	IteratingLogProcess.cpp
	$(ECELL3_DMC) -o IteratingLogProcess.so --ldflags=SpatiocyteProcess.so IteratingLogProcess.cpp

CoordinateLogProcess.so: 	CoordinateLogProcess.cpp
	$(ECELL3_DMC) -o CoordinateLogProcess.so --ldflags="SpatiocyteProcess.so IteratingLogProcess.so" CoordinateLogProcess.cpp

DiffusionInfluencedReactionProcess.so: 	DiffusionInfluencedReactionProcess.cpp
	$(ECELL3_DMC) -o DiffusionInfluencedReactionProcess.so --ldflags=ReactionProcess.so DiffusionInfluencedReactionProcess.cpp

OscillationAnalysisProcess.so: 	OscillationAnalysisProcess.cpp
	$(ECELL3_DMC) -o OscillationAnalysisProcess.so --ldflags="SpatiocyteProcess.so" OscillationAnalysisProcess.cpp

SpatiocyteNextReactionProcess.so: 	SpatiocyteNextReactionProcess.cpp
	$(ECELL3_DMC) -o SpatiocyteNextReactionProcess.so --ldflags=ReactionProcess.so SpatiocyteNextReactionProcess.cpp

MoleculePopulateProcess.so: 	MoleculePopulateProcess.cpp
	$(ECELL3_DMC) -o MoleculePopulateProcess.so --ldflags="SpatiocyteProcess.so DiffusionInfluencedReactionProcess.so" MoleculePopulateProcess.cpp

DiffusionProcess.so: 	DiffusionProcess.cpp
	$(ECELL3_DMC) -o DiffusionProcess.so --ldflags="SpatiocyteProcess.so DiffusionInfluencedReactionProcess.so" DiffusionProcess.cpp

PolymerizationParameterProcess.so: 	PolymerizationParameterProcess.cpp
	$(ECELL3_DMC) -o PolymerizationParameterProcess.so --ldflags="SpatiocyteProcess.so" PolymerizationParameterProcess.cpp

PolymerizationProcess.so: 	PolymerizationProcess.cpp
	$(ECELL3_DMC) -o PolymerizationProcess.so --ldflags="SpatiocyteProcess.so DiffusionInfluencedReactionProcess.so" PolymerizationProcess.cpp

PeriodicBoundaryDiffusionProcess.so: 	PeriodicBoundaryDiffusionProcess.cpp
	$(ECELL3_DMC) -o PeriodicBoundaryDiffusionProcess.so --ldflags=DiffusionProcess.so PeriodicBoundaryDiffusionProcess.cpp

PolymerFragmentationProcess.so: 	PolymerFragmentationProcess.cpp
	$(ECELL3_DMC) -o PolymerFragmentationProcess.so --ldflags=ReactionProcess.so PolymerFragmentationProcess.cpp



%.so: %.cpp
	$(ECELL3_DMC) $< 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

$(SPATIOCYTE):$(OBJECTS)
	$(CXX) -v -o $@ $(OBJECTS) $(GUILIBS)

gui:	$(SPATIOCYTE)

back:
	mkdir ../backup/$(VER)
	cp *.?pp ../backup/$(VER)
	cp Makefile ../backup/$(VER)
	cp *.em ../backup/$(VER)
	cp *.txt ../backup/$(VER)
	cp *.ess ../backup/$(VER)
	cp *.py ../backup/$(VER)
#	cp *.csv ../backup/$(VER)

backBench:
	mkdir ../backup/$(VER)
	cp *.?pp ../backup/$(VER)
	cp Makefile ../backup/$(VER)
	cp *.em ../backup/$(VER)
	cp *.ess ../backup/$(VER)
clean: 
	rm -f *.so

allclean: 
	rm -f *.so *.dat *.o *.png $(SPATIOCYTE)

pclean: 
	rm -f *Process.so

p:
	rm -f *Process.so
	make -j
	gecell -f membrane.eml
