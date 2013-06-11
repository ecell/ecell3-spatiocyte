import math

Filaments = 13
RotateAngle = math.pi
MTRadius = 12.5e-9
VoxelRadius = 0.8e-8
KinesinRadius = 0.4e-8
neuriteRadius = 0.15e-6
neuriteLength = 5e-6
somaRadius = 1.0e-6
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
react.k = 2.5863133e-22
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

def rotatePointAlongVector(P, C, N, angle):
  x = P[0]
  y = P[1]
  z = P[2]
  a = C[0]
  b = C[1]
  c = C[2]
  u = N[0]
  v = N[1]
  w = N[2]
  u2 = u*u
  v2 = v*v
  w2 = w*w
  cosT = math.cos(angle)
  oneMinusCosT = 1-cosT
  sinT = math.sin(angle)
  xx = (a*(v2 + w2) - u*(b*v + c*w - u*x - v*y - w*z)) * oneMinusCosT + x*cosT + (-c*v + b*w - w*y + v*z)*sinT
  yy = (b*(u2 + w2) - v*(a*u + c*w - u*x - v*y - w*z)) * oneMinusCosT + y*cosT + (c*u - a*w + w*x - u*z)*sinT
  zz = (c*(u2 + v2) - w*(a*u + b*v - u*x - v*y - w*z)) * oneMinusCosT + z*cosT + (-b*u + a*v - v*x + u*y)*sinT
  return [xx, yy, zz]

somaMTs = 16
somaMTrotateAngle = math.pi*2/somaMTs
somaMTorigin = [0.5, 0.0, 0.0]
somaMTvectorZ = [0.0, 0.0, 1.0]
somaMTvectorZpoint = [0.0, 0.0, 0.0]
for i in range(somaMTs):
  Microtubule = theSimulator.createEntity('MicrotubuleProcess', 'Process:/Soma:Microtubule%d' %(i))
  P = rotatePointAlongVector(somaMTorigin, somaMTvectorZpoint, somaMTvectorZ, somaMTrotateAngle*i)
  Microtubule.OriginX = P[0]
  Microtubule.OriginY = P[1]
  Microtubule.OriginZ = P[2]
  Microtubule.RotateX = 0
  Microtubule.RotateY = 0
  Microtubule.RotateZ =  somaMTrotateAngle*i
  Microtubule.Radius = MTRadius
  Microtubule.SubunitRadius = KinesinRadius
  Microtubule.Length = somaRadius*0.7
  Microtubule.Filaments = Filaments
  Microtubule.Periodic = 0
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesin' ]]
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:MTKinesinATP' ]]
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:actTubulin' ]]
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:Tubulin' , '-1']]
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:TubulinM' , '-2']]
  Microtubule.VariableReferenceList = [['_', 'Variable:/Soma:TubulinP' , '-3']]

neuritesLengthX = [neuriteLength, neuriteLength, neuriteLength, neuriteLength]
neuritesOriginX = [0.55, -0.55]
neuritesOriginY = [0, 0]
neuritesRotateZ = [0, math.pi]

MTsOriginY = [0, 0.6, -0.6, 0, 0]
MTsOriginZ = [0, 0, 0, 0.6, -0.6]

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

  for j in range(5):
    Microtubule = theSimulator.createEntity('MicrotubuleProcess', 'Process:/neurite%d:Microtubule%d' %(i, j))
    Microtubule.OriginX = 0
    Microtubule.OriginY = MTsOriginY[j]
    Microtubule.OriginZ = MTsOriginZ[j]
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

run(10)

