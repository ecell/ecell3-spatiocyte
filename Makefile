OBJS=\
SpatiocyteVisualizer

DMS=\
SpatiocyteStepper\
SpatiocyteProcess\
MassActionProcess\
ReactionProcess\
VisualizationLogProcess\
MicrotubuleProcess\
HistogramLogProcess\
H5VisualizationLogProcess\
IteratingLogProcess\
DiffusionInfluencedReactionProcess\
OscillationAnalysisProcess\
SpatiocyteNextReactionProcess\
MicroscopyTrackingProcess\
CoordinateLogProcess\
MoleculePopulateProcess\
DiffusionProcess\
PolymerFragmentationProcess\
PolymerizationParameterProcess\
PolymerizationProcess\
PeriodicBoundaryDiffusionProcess

ECELL3_DMC = ecell3-dmc
CXX = g++
CXXFLAGS = -Wall -O3 -g
CXXFLAGS += $(shell pkg-config --cflags gtkmm-2.4 gtkglextmm-x11-1.2)
CPPFLAGS = -DG_DISABLE_DEPRECATED -DGDK_PIXBUF_DISABLE_DEPRECATED -DPNG_SKIP_SETJMP_CHECK # -DGDK_DISABLE_DEPRECATED 
GUILIBS =
GUILIBS += $(shell pkg-config --libs gtkmm-2.4 gtkglextmm-x11-1.2 libpng)
SPATIOCYTE = spatiocyte
OBJECTS=${OBJS:=.o}
SOS=${DMS:=.so}

all:	$(SOS) $(SPATIOCYTE)

VisualizationLogProcess.so: 	VisualizationLogProcess.cpp
	$(ECELL3_DMC) -o VisualizationLogProcess.so --ldflags=SpatiocyteProcess.so VisualizationLogProcess.cpp

H5VisualizationLogProcess.so: 	H5VisualizationLogProcess.cpp
	$(ECELL3_DMC) -o H5VisualizationLogProcess.so --ldflags=SpatiocyteProcess.so --ldflags=-lhdf5 --ldflags=-lhdf5_cpp H5VisualizationLogProcess.cpp

MicroscopyTrackingProcess.so: 	MicroscopyTrackingProcess.cpp
	$(ECELL3_DMC) -o MicroscopyTrackingProcess.so --ldflags="SpatiocyteProcess.so VisualizationLogProcess.so" MicroscopyTrackingProcess.cpp

ReactionProcess.so: 	ReactionProcess.cpp
	$(ECELL3_DMC) -o ReactionProcess.so --ldflags=SpatiocyteProcess.so ReactionProcess.cpp

IteratingLogProcess.so: 	IteratingLogProcess.cpp
	$(ECELL3_DMC) -o IteratingLogProcess.so --ldflags=SpatiocyteProcess.so IteratingLogProcess.cpp

CoordinateLogProcess.so: 	CoordinateLogProcess.cpp
	$(ECELL3_DMC) -o CoordinateLogProcess.so --ldflags="SpatiocyteProcess.so IteratingLogProcess.so" CoordinateLogProcess.cpp

DiffusionInfluencedReactionProcess.so: 	DiffusionInfluencedReactionProcess.cpp
	$(ECELL3_DMC) -o DiffusionInfluencedReactionProcess.so --ldflags=ReactionProcess.so DiffusionInfluencedReactionProcess.cpp

MicrotubuleProcess.so: 	MicrotubuleProcess.cpp
	$(ECELL3_DMC) -o MicrotubuleProcess.so --ldflags="SpatiocyteProcess.so" MicrotubuleProcess.cpp

HistogramLogProcess.so: 	HistogramLogProcess.cpp
	$(ECELL3_DMC) -o HistogramLogProcess.so --ldflags="SpatiocyteProcess.so" HistogramLogProcess.cpp

OscillationAnalysisProcess.so: 	OscillationAnalysisProcess.cpp
	$(ECELL3_DMC) -o OscillationAnalysisProcess.so --ldflags="SpatiocyteProcess.so" OscillationAnalysisProcess.cpp

SpatiocyteNextReactionProcess.so: 	SpatiocyteNextReactionProcess.cpp
	$(ECELL3_DMC) -o SpatiocyteNextReactionProcess.so --ldflags=ReactionProcess.so SpatiocyteNextReactionProcess.cpp

MoleculePopulateProcess.so: 	MoleculePopulateProcess.cpp
	$(ECELL3_DMC) -o MoleculePopulateProcess.so --ldflags="SpatiocyteProcess.so DiffusionInfluencedReactionProcess.so" MoleculePopulateProcess.cpp

DiffusionProcess.so: 	DiffusionProcess.cpp
	$(ECELL3_DMC) -o DiffusionProcess.so --ldflags="SpatiocyteProcess.so DiffusionInfluencedReactionProcess.so" DiffusionProcess.cpp

CompartmentGrowthProcess.so: 	CompartmentGrowthProcess.cpp
	$(ECELL3_DMC) -o CompartmentGrowthProcess.so --ldflags="SpatiocyteProcess.so" CompartmentGrowthProcess.cpp

PolymerizationParameterProcess.so: 	PolymerizationParameterProcess.cpp
	$(ECELL3_DMC) -o PolymerizationParameterProcess.so --ldflags="SpatiocyteProcess.so" PolymerizationParameterProcess.cpp

PolymerizationProcess.so: 	PolymerizationProcess.cpp
	$(ECELL3_DMC) -o PolymerizationProcess.so --ldflags="SpatiocyteProcess.so DiffusionInfluencedReactionProcess.so PolymerFragmentationProcess.so" PolymerizationProcess.cpp

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

clean: 
	rm -f *.so *.o $(SPATIOCYTE)
