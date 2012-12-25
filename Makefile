OBJS=\
SpatiocyteVisualizer

DMS=\
SpatiocyteStepper\
SpatiocyteProcess\
VisualizationLogProcess\
ReactionProcess\
DiffusionInfluencedReactionProcess\
MoleculePopulateProcess\
DiffusionProcess
#\
DMS=\
SpatiocyteStepper\
SpatiocyteProcess\
MassActionProcess\
ReactionProcess\
VisualizationLogProcess\
CompartmentProcess\
ErythrocyteProcess\
H5VisualizationLogProcess\
DiffusionInfluencedReactionProcess\
IteratingLogProcess\
OscillationAnalysisProcess\
SpatiocyteNextReactionProcess\
MicroscopyTrackingProcess\
DiffusionProcess\
MoleculeLogProcess\
HistogramLogProcess\
MoleculePopulateProcess\
TagProcess\
PeriodicBoundaryDiffusionProcess
#PolymerFragmentationProcess\
MicrotubuleProcess\
PolymerizationParameterProcess\
PolymerizationProcess

CXXFLAGS = -Wall -O3 -g
ECELL3_DMC = ecell3-dmc #--no-define-cxxflags="$(CXXFLAGS)"
CXX = g++
GUIFLAGS = $(CXXFLAGS) $(shell pkg-config --cflags gtkmm-2.4 gtkglextmm-x11-1.2)
CPPFLAGS = -DG_DISABLE_DEPRECATED -DGDK_PIXBUF_DISABLE_DEPRECATED -DPNG_SKIP_SETJMP_CHECK # -DGDK_DISABLE_DEPRECATED 
GUILIBS =
GUILIBS += $(shell pkg-config --libs gtkmm-2.4 gtkglextmm-x11-1.2 libpng)
SPATIOCYTE = spatiocyte
OBJECTS=${OBJS:=.o}
SOS=${DMS:=.so}

all:	$(SPATIOCYTE) $(SOS) 

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

MoleculeLogProcess.so: 	MoleculeLogProcess.cpp
	$(ECELL3_DMC) -o MoleculeLogProcess.so --ldflags="SpatiocyteProcess.so IteratingLogProcess.so" MoleculeLogProcess.cpp

DiffusionInfluencedReactionProcess.so: 	DiffusionInfluencedReactionProcess.cpp
	$(ECELL3_DMC) -o DiffusionInfluencedReactionProcess.so --ldflags=ReactionProcess.so DiffusionInfluencedReactionProcess.cpp

MicrotubuleProcess.so: 	MicrotubuleProcess.cpp
	$(ECELL3_DMC) -o MicrotubuleProcess.so --ldflags="SpatiocyteProcess.so" MicrotubuleProcess.cpp

CompartmentProcess.so: 	CompartmentProcess.cpp
	$(ECELL3_DMC) -o CompartmentProcess.so --ldflags="SpatiocyteProcess.so" CompartmentProcess.cpp

ErythrocyteProcess.so: 	ErythrocyteProcess.cpp
	$(ECELL3_DMC) -o ErythrocyteProcess.so --ldflags="SpatiocyteProcess.so" ErythrocyteProcess.cpp

HistogramLogProcess.so: 	HistogramLogProcess.cpp
	$(ECELL3_DMC) -o HistogramLogProcess.so --ldflags="SpatiocyteProcess.so IteratingLogProcess.so" HistogramLogProcess.cpp

OscillationAnalysisProcess.so: 	OscillationAnalysisProcess.cpp
	$(ECELL3_DMC) -o OscillationAnalysisProcess.so --ldflags="SpatiocyteProcess.so" OscillationAnalysisProcess.cpp

SpatiocyteNextReactionProcess.so: 	SpatiocyteNextReactionProcess.cpp
	$(ECELL3_DMC) -o SpatiocyteNextReactionProcess.so --ldflags=ReactionProcess.so SpatiocyteNextReactionProcess.cpp

TagProcess.so: 	TagProcess.cpp
	$(ECELL3_DMC) -o TagProcess.so --ldflags="SpatiocyteProcess.so" TagProcess.cpp

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

dms: $(SOS) 

%.so: %.cpp
	$(ECELL3_DMC) $< 

%.o: %.cpp
	$(CXX) $(GUIFLAGS) $(CPPFLAGS) -c -o $@ $<

$(SPATIOCYTE):$(OBJECTS)
	$(CXX) -v -o $@ $(OBJECTS) $(GUILIBS)

gui:	$(SPATIOCYTE)

clean: 
	rm -f *.so *.o $(SPATIOCYTE)
