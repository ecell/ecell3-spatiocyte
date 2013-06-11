import math

Filaments = 13
RotateAngle = math.pi
MTRadius = 12.5e-9
VoxelRadius = 0.8e-8
KinesinRadius = 0.4e-8
neuriteRadius = 0.15e-6
neuriteLength = 5e-6
somaRadius = 0.75e-6
MTLength = neuriteLength*0.9

singleMTVolumeVoxels = 717256.0
singleNeuriteVolumeVoxels = 48325789.0
totalKinesins = 1000*singleMTVolumeVoxels/singleNeuriteVolumeVoxels
print "totalKinesins:", totalKinesins

theSimulator.createStepper('SpatiocyteStepper', 'SS').VoxelRadius = VoxelRadius

# Create the root container compartment using the default Cuboid geometry:
theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = somaRadius*2+neuriteLength*2
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = somaRadius*2.1
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = neuriteRadius*4.1
theSimulator.createEntity('Variable', 'Variable:/:VACANT')

# Create the Soma compartment of the Neuron:
theSimulator.createEntity('System', 'System:/:Soma').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Soma:GEOMETRY').Value = 1
theSimulator.createEntity('Variable', 'Variable:/Soma:LENGTHX').Value = somaRadius*2
theSimulator.createEntity('Variable', 'Variable:/Soma:LENGTHY').Value = somaRadius*2
theSimulator.createEntity('Variable', 'Variable:/Soma:LENGTHZ').Value = neuriteRadius*4
theSimulator.createEntity('Variable', 'Variable:/Soma:VACANT')
theSimulator.createEntity('Variable', 'Variable:/Soma:Kinesin').Value = 60
theSimulator.createEntity('Variable', 'Variable:/Soma:MTKinesin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/Soma:MTKinesinATP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/Soma:Tubulin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/Soma:actTubulin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/Soma:TubulinM' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/Soma:TubulinP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/Soma:GFP' ).Value = 0
# Create the Soma membrane:
theSimulator.createEntity('System', 'System:/Soma:Membrane').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Soma/Membrane:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Soma/Membrane:VACANT')

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/Soma:populate')
populator.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin']]
populator.VariableReferenceList = [['_', 'Variable:/Soma:actTubulin']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/Soma:populateK')
populator.VariableReferenceList = [['_', 'Variable:/Soma:Kinesin']]
populator.OriginX = -1.0
populator.UniformRadiusX = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/Soma:detachPlus')
react.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:TubulinP','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:Kinesin','1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:TubulinP','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/Soma:explicitAttach')
react.VariableReferenceList = [['_', 'Variable:/Soma:Kinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:Tubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin','1']]
react.k = 2.5863133e-24
#react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/Soma:explicitAttachAct')
react.VariableReferenceList = [['_', 'Variable:/Soma:Kinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:actTubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin','1']]
react.k = 1.5863133e-21
#react.p = 1

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/Soma:diffuseKinesin')
diffuse.VariableReferenceList = [['_', 'Variable:/Soma:Kinesin']]
diffuse.D = 4e-12

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/Soma:detach')
react.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:actTubulin','1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:Kinesin','1']]
react.SearchVacant = 1
react.k = 0.25

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/Soma:hydrolysis')
react.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin','1']]
react.SearchVacant = 1
react.k = 100

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/Soma:phosphorylate')
react.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesinATP','1']]
react.SearchVacant = 1
react.k = 145

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/Soma:ratchet')
react.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesinATP','1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:actTubulin','1']]
react.BindingSite = 1
react.k = 55

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/Soma:inactive')
react.VariableReferenceList = [['_', 'Variable:/Soma:actTubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/Soma:Tubulin','1']]
react.k = 5

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/Soma:diffusePlus')
diffuse.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin']]
diffuse.VariableReferenceList = [['_', 'Variable:/Soma:actTubulin', '1']]
diffuse.D = 0.04e-12

#tagger = theSimulator.createEntity('TagProcess', 'Process:/Soma:tagger')
#tagger.VariableReferenceList = [['_', 'Variable:/Soma:GFP', '-1' ]]
#tagger.VariableReferenceList = [['_', 'Variable:/Soma:Kinesin', '5' ]]
#tagger.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin']]
#tagger.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesinATP']]

#life = theSimulator.createEntity('LifetimeLogProcess', 'Process:/Soma:lifetime')
#life.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin', '-1']]
#life.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesinATP', '-1']]
#life.VariableReferenceList = [['_', 'Variable:/Soma:Kinesin', '1']]

visualLogger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/Soma:visualLogger')
#visualLogger.VariableReferenceList = [['_', 'Variable:/Soma/Membrane:VACANT']]
visualLogger.VariableReferenceList = [['_', 'Variable:/Soma:Tubulin']]
visualLogger.VariableReferenceList = [['_', 'Variable:/Soma:TubulinM']]
visualLogger.VariableReferenceList = [['_', 'Variable:/Soma:TubulinP']]
visualLogger.VariableReferenceList = [['_', 'Variable:/Soma:Kinesin', '10600']]
visualLogger.VariableReferenceList = [['_', 'Variable:/Soma:actTubulin']]
visualLogger.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin', '10600']]
visualLogger.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesinATP', '10600']]
visualLogger.LogInterval = 1

neuritesLengthX = [neuriteLength, neuriteLength, neuriteLength, neuriteLength]
neuritesOriginX = [0.5, -0.5] 
neuritesOriginY = [0, 0]
neuritesRotateZ = [0, math.pi]
for i in range(2):
  # Create the neurite:
  theSimulator.createEntity('System', 'System:/:neurite%d' %i).StepperID = 'SS'
  theSimulator.createEntity('Variable', 'Variable:/neurite%d:GEOMETRY' %i).Value = 3
  theSimulator.createEntity('Variable', 'Variable:/neurite%d:LENGTHX' %i).Value = neuritesLengthX[i]
  theSimulator.createEntity('Variable', 'Variable:/neurite%d:LENGTHY' %i).Value = neuriteRadius*2
  theSimulator.createEntity('Variable', 'Variable:/neurite%d:ORIGINX' %i).Value = neuritesOriginX[i]
  theSimulator.createEntity('Variable', 'Variable:/neurite%d:ORIGINY' %i).Value = neuritesOriginY[i]
  theSimulator.createEntity('Variable', 'Variable:/neurite%d:ORIGINZ' %i).Value = 0
  theSimulator.createEntity('Variable', 'Variable:/neurite%d:ROTATEZ' %i).Value = neuritesRotateZ[i]
  theSimulator.createEntity('Variable', 'Variable:/neurite%d:VACANT' %i)
  theSimulator.createEntity('Variable', 'Variable:/neurite%d:DIFFUSIVE' %i).Name = '/:Soma'
  # Create the neurite membrane:
  theSimulator.createEntity('System', 'System:/neurite%d:Membrane' %i).StepperID = 'SS'
  theSimulator.createEntity('Variable', 'Variable:/neurite%d/Membrane:DIMENSION' %i).Value = 2
  theSimulator.createEntity('Variable', 'Variable:/neurite%d/Membrane:VACANT' %i)
  theSimulator.createEntity('Variable', 'Variable:/neurite%d/Membrane:DIFFUSIVE' %i).Name = '/Soma:Membrane'

  Microtubule = theSimulator.createEntity('MicrotubuleProcess', 'Process:/neurite%d:Microtubule' %i)
  Microtubule.OriginX = 0
  Microtubule.OriginY = 0
  Microtubule.OriginZ = 0
  Microtubule.RotateX = 0
  Microtubule.RotateY = 0
  Microtubule.RotateZ =  0
  Microtubule.Radius = MTRadius
  Microtubule.SubunitRadius = KinesinRadius
  Microtubule.Length = MTLength
  Microtubule.Filaments = Filaments
  Microtubule.Periodic = 0
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin' ]]
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesinATP' ]]
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:actTubulin' ]]
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:Tubulin' , '-1']]
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:TubulinM' , '-2']]
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:TubulinP' , '-3']]

#Microtubule = theSimulator.createEntity('MicrotubuleProcess', 'Process:/Soma:MicrotubuleLeft')
#Microtubule.OriginX = -0.5
#Microtubule.OriginY = 0
#Microtubule.OriginZ = 0
#Microtubule.RotateX = 0
#Microtubule.RotateY = 0
#Microtubule.RotateZ = RotateAngle
#Microtubule.Radius = MTRadius
#Microtubule.SubunitRadius = KinesinRadius
#Microtubule.Length = MTLength
#Microtubule.Filaments = Filaments
#Microtubule.Periodic = 0
#Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin' ]]
#Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesinATP' ]]
#Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:actTubulin' ]]
#Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:Tubulin' , '-1']]
#Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:TubulinM' , '-2']]
#Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:TubulinP' , '-3']]

run(2000)

