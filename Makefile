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
FluorescentProteinImagingProcess\
CoordinateLogProcess\
MoleculePopulateProcess\
DiffusionProcess\
PolymerizationParameterProcess\
PolymerizationProcess\
PeriodicBoundaryDiffusionProcess\
PolymerFragmentationProcess

CXX = g++
CXXFLAGS = -Wall -O3 -I/usr/include/freetype2 -I/usr/include/gtkglextmm-1.2 -I/usr/lib/gtkglextmm-1.2/include -I/usr/include/gtkglext-1.0 -I/usr/lib/gtkglext-1.0/include -I/usr/include/gtk-2.0 -I/usr/lib/gtk-2.0/include -I/usr/include/pango-1.0 -I/usr/include/cairo -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include -I/usr/include/gdkmm-2.4 -I/usr/lib/gdkmm-2.4/include -I/usr/include/glibmm-2.4 -I/usr/include/giomm-2.4 -I/usr/lib/glibmm-2.4/include -I/usr/include/cairomm-1.0 -I/usr/include/pangomm-1.4 -I/usr/lib/pangomm-1.4/include -I/usr/include/sigc++-2.0 -I/usr/lib/sigc++-2.0/include -I/usr/include/atk-1.0 -I/usr/include/gtkmm-2.4 -I/usr/lib/gtkmm-2.4/include -I/usr/include/atkmm-1.6 -DG_DISABLE_DEPRECATED -DGDK_DISABLE_DEPRECATED -DGDK_PIXBUF_DISABLE_DEPRECATED -DPNG_SKIP_SETJMP_CHECK -g
GUILIBS = -lgtkglextmm-x11-1.2 -lgtkglext-x11-1.0 -lgdkglext-x11-1.0 -lGLU -lGL -lXmu -lXt -lSM -lICE -lpangox-1.0 -lgtkmm-2.4 -lgdkmm-2.4 -lgtk-x11-2.0 -lpangomm-1.4 -lglibmm-2.4 -lsigc-2.0 -lgdk-x11-2.0 -latk-1.0 -lgdk_pixbuf-2.0 -lm -lpangocairo-1.0 -lpango-1.0 -lcairo -lgobject-2.0 -lgmodule-2.0 -ldl -lglib-2.0 
SPATIOCYTE = spatiocyte
OBJECTS=${OBJS:=.o}
SOS=${DMS:=.so}

all:	$(SOS) $(SPATIOCYTE)

VisualizationLogProcess.so: 	VisualizationLogProcess.cpp
	ecell3-dmc -o VisualizationLogProcess.so --ldflags=SpatiocyteProcess.so VisualizationLogProcess.cpp

FluorescentProteinImagingProcess.so: 	FluorescentProteinImagingProcess.cpp
	ecell3-dmc -o FluorescentProteinImagingProcess.so --ldflags="SpatiocyteProcess.so VisualizationLogProcess.so" FluorescentProteinImagingProcess.cpp

ReactionProcess.so: 	ReactionProcess.cpp
	ecell3-dmc -o ReactionProcess.so --ldflags=SpatiocyteProcess.so ReactionProcess.cpp

IteratingLogProcess.so: 	IteratingLogProcess.cpp
	ecell3-dmc -o IteratingLogProcess.so --ldflags=SpatiocyteProcess.so IteratingLogProcess.cpp

CoordinateLogProcess.so: 	CoordinateLogProcess.cpp
	ecell3-dmc -o CoordinateLogProcess.so --ldflags="SpatiocyteProcess.so IteratingLogProcess.so" CoordinateLogProcess.cpp

DiffusionInfluencedReactionProcess.so: 	DiffusionInfluencedReactionProcess.cpp
	ecell3-dmc -o DiffusionInfluencedReactionProcess.so --ldflags=ReactionProcess.so DiffusionInfluencedReactionProcess.cpp

OscillationAnalysisProcess.so: 	OscillationAnalysisProcess.cpp
	ecell3-dmc -o OscillationAnalysisProcess.so --ldflags="SpatiocyteProcess.so" OscillationAnalysisProcess.cpp

SpatiocyteNextReactionProcess.so: 	SpatiocyteNextReactionProcess.cpp
	ecell3-dmc -o SpatiocyteNextReactionProcess.so --ldflags=ReactionProcess.so SpatiocyteNextReactionProcess.cpp

MoleculePopulateProcess.so: 	MoleculePopulateProcess.cpp
	ecell3-dmc -o MoleculePopulateProcess.so --ldflags="SpatiocyteProcess.so DiffusionInfluencedReactionProcess.so" MoleculePopulateProcess.cpp

DiffusionProcess.so: 	DiffusionProcess.cpp
	ecell3-dmc -o DiffusionProcess.so --ldflags="SpatiocyteProcess.so DiffusionInfluencedReactionProcess.so" DiffusionProcess.cpp

PolymerizationParameterProcess.so: 	PolymerizationParameterProcess.cpp
	ecell3-dmc -o PolymerizationParameterProcess.so --ldflags="SpatiocyteProcess.so" PolymerizationParameterProcess.cpp

PolymerizationProcess.so: 	PolymerizationProcess.cpp
	ecell3-dmc -o PolymerizationProcess.so --ldflags="SpatiocyteProcess.so DiffusionInfluencedReactionProcess.so" PolymerizationProcess.cpp

PeriodicBoundaryDiffusionProcess.so: 	PeriodicBoundaryDiffusionProcess.cpp
	ecell3-dmc -o PeriodicBoundaryDiffusionProcess.so --ldflags=DiffusionProcess.so PeriodicBoundaryDiffusionProcess.cpp

PolymerFragmentationProcess.so: 	PolymerFragmentationProcess.cpp
	ecell3-dmc -o PolymerFragmentationProcess.so --ldflags=ReactionProcess.so PolymerFragmentationProcess.cpp



%.so: 	%.cpp
	ecell3-dmc $< 

$(SPATIOCYTE):	$(OBJECTS)
	$(CXX) -v -o spatiocyte $(OBJECTS) $(GUILIBS)

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
