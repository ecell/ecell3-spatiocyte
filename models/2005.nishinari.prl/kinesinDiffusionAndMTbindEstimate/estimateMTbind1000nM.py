import math
import random
minDist = 75e-9
dendriteRadius = 0.5e-6
Protofilaments = 13
RotateAngle = 0 #math.pi/4
MTRadius = 12.5e-9
VoxelRadius = 0.4e-8
KinesinRadius = 0.4e-8

mtOriginX = [0]#, 0]
mtOriginZ = [0]#, 0.5]
mtOriginY = [0]#, 0.5]
completedLengths = [KinesinRadius*2*700]#, 0.5e-6]
dendriteLength = completedLengths[0]*1.2

theSimulator.createStepper('SpatiocyteStepper', 'SS').VoxelRadius = VoxelRadius
theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 3
theSimulator.createEntity('Variable', 'Variable:/:ROTATEZ').Value = RotateAngle
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = dendriteLength
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = dendriteRadius*2
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:Kinesin').Value = 3071
theSimulator.createEntity('Variable', 'Variable:/:MTKinesin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MTKinesinATP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Tubulin' ).Value = 0

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populate')
populator.VariableReferenceList = [['_', 'Variable:/:Kinesin']]

visualLogger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:visualLogger')
visualLogger.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:Tubulin' ]]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
visualLogger.LogInterval = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttach')
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','0']]
react.p = 1

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseKinesin')
diffuse.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
diffuse.D = 8e-12

for i in range(len(completedLengths)):
  Microtubule = theSimulator.createEntity('MicrotubuleProcess', 'Process:/:Microtubule%d' %i)
  Microtubule.OriginX = mtOriginX[i]
  Microtubule.OriginY = mtOriginY[i]
  Microtubule.OriginZ = mtOriginZ[i]
  Microtubule.RotateX = 0
  Microtubule.RotateY = 0
  Microtubule.RotateZ = RotateAngle
  Microtubule.Length = completedLengths[i]
  Microtubule.Radius = MTRadius
  Microtubule.SubunitRadius = KinesinRadius
  Microtubule.Filaments = Protofilaments
  Microtubule.Periodic = 1
  Microtubule.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
  Microtubule.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
  Microtubule.VariableReferenceList = [['_', 'Variable:/:Tubulin' , '-1']]
  Iterator = theSimulator.createEntity('IteratingLogProcess', 'Process:/:Iterator%d' %i)
  Iterator.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
  Iterator.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
  Iterator.VariableReferenceList = [['_', 'Variable:/:Kinesin' ]]
  Iterator.LogEnd = 0.5
  Iterator.LogInterval = 0.01
  Iterator.Iterations = 1
  Iterator.FileName = "estimate1000nM.csv"
run(Iterator.LogEnd + 0.1)


