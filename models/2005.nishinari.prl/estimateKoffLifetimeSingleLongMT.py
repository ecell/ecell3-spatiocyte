import math

Filaments = 13
RotateAngle = 0 #math.pi/2
MTRadius = 12.5e-9
VoxelRadius = 0.4e-8
KinesinRadius = 0.4e-8
dendriteRadius = 0.35e-6
dendriteLength = 60e-6
totalMTLength = 8e-6

theSimulator.createStepper('SpatiocyteStepper', 'SS').VoxelRadius = VoxelRadius
theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 2
theSimulator.createEntity('Variable', 'Variable:/:ROTATEZ').Value = RotateAngle
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = dendriteLength
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = dendriteRadius*2
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:Kinesin').Value = 10 
theSimulator.createEntity('Variable', 'Variable:/:MTKinesin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MTKinesinATP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Tubulin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:TubulinM' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:TubulinP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:GFP' ).Value = 0

theSimulator.createEntity('System', 'System:/:Membrane').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Membrane:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Membrane:VACANT')

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populate')
populator.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populateK')
populator.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
populator.OriginX = -0.7

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttach')
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:detachPlus')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinP','1']]
react.p = 1

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseKinesin')
diffuse.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
diffuse.D = 1.8e-12

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:detach')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','1']]
react.SearchVacant = 1
react.k = 0.25

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:hydrolysis')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.SearchVacant = 1
react.k = 100

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:phosphorylate')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','1']]
react.SearchVacant = 1
react.k = 145

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:ratchet')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','1']]
react.BindingSite = 1
react.k = 55

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePlus')
diffuse.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
diffuse.D = 0.04e-12

#tagger = theSimulator.createEntity('TagProcess', 'Process:/:tagger')
#tagger.VariableReferenceList = [['_', 'Variable:/:GFP', '-1' ]]
#tagger.VariableReferenceList = [['_', 'Variable:/:Kinesin', '10' ]]
#tagger.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
#tagger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP']]

life = theSimulator.createEntity('LifetimeLogProcess', 'Process:/:lifetime')
life.VariableReferenceList = [['_', 'Variable:/:MTKinesin', '-1']]
life.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP', '-1']]
life.VariableReferenceList = [['_', 'Variable:/:Kinesin', '1']]
life.Iterations = 500
life.LogEnd = 99

#micro = theSimulator.createEntity('MicroscopyTrackingProcess', 'Process:/:track')
#micro.VariableReferenceList = [['_', 'Variable:/:MTKinesin', '1']]
#micro.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP', '1']]
#micro.VariableReferenceList = [['_', 'Variable:/:MTKinesin', '-1']]
#micro.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP', '-1']]
#micro.ExposureTime = 0.03

#visualLogger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:visualLogger')
#visualLogger.VariableReferenceList = [['_', 'Variable:/:Tubulin']]
#visualLogger.VariableReferenceList = [['_', 'Variable:/:TubulinM']]
#visualLogger.VariableReferenceList = [['_', 'Variable:/:TubulinP']]
#visualLogger.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
#visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
#visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
#visualLogger.VariableReferenceList = [['_', 'Variable:/:GFP' ]]
#visualLogger.LogInterval = 1e-2

Microtubule = theSimulator.createEntity('MicrotubuleProcess', 'Process:/:Microtubule')
Microtubule.OriginX = 0
Microtubule.OriginY = 0
Microtubule.OriginZ = 0
Microtubule.RotateX = 0
Microtubule.RotateY = 0
Microtubule.RotateZ = RotateAngle
Microtubule.Radius = MTRadius
Microtubule.SubunitRadius = KinesinRadius
Microtubule.Length = 58e-6
Microtubule.Filaments = 13
Microtubule.Periodic = 0
Microtubule.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
Microtubule.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
Microtubule.VariableReferenceList = [['_', 'Variable:/:Tubulin' , '-1']]
Microtubule.VariableReferenceList = [['_', 'Variable:/:TubulinM' , '-2']]
Microtubule.VariableReferenceList = [['_', 'Variable:/:TubulinP' , '-3']]

run(100)

