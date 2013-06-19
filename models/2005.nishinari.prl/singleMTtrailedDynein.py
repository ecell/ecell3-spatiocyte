Filaments = 13
RotateAngle = 0 #math.pi/4
MTRadius = 12.5e-9
VoxelRadius = 0.4e-8
KinesinRadius = 0.4e-8
dendriteRadius = 0.15e-6
dendriteLength = 1.4e-6
MTLength = dendriteLength*0.95

theSimulator.createStepper('SpatiocyteStepper', 'SS').VoxelRadius = VoxelRadius
theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 2
theSimulator.createEntity('Variable', 'Variable:/:ROTATEZ').Value = RotateAngle
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = dendriteLength
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = dendriteRadius*2
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:Kinesin').Value = 10
theSimulator.createEntity('Variable', 'Variable:/:Dynein').Value = 10
theSimulator.createEntity('Variable', 'Variable:/:MTKinesin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MTDynein' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MTKinesinATP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MTDyneinATP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Tubulin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:actTubulin' ).Value = 2000
theSimulator.createEntity('Variable', 'Variable:/:TubulinM' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:TubulinP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:GFP' ).Value = 0

theSimulator.createEntity('System', 'System:/:Membrane').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Membrane:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Membrane:VACANT')

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populate')
populator.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
populator.VariableReferenceList = [['_', 'Variable:/:MTDynein']]
populator.VariableReferenceList = [['_', 'Variable:/:actTubulin']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populateK')
populator.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
populator.VariableReferenceList = [['_', 'Variable:/:Dynein']]
populator.OriginZ = -1.5
populator.OriginY = 1.5

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:detachPlus')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinP','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:detachMinus')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinM','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinM','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:detachPlusDynein')
react.VariableReferenceList = [['_', 'Variable:/:MTDynein','-1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Dynein','1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinP','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:detachMinusDynein')
react.VariableReferenceList = [['_', 'Variable:/:MTDynein','-1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinM','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Dynein','1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinM','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttach')
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.p = 0.001

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttachAct')
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttachDynein')
react.VariableReferenceList = [['_', 'Variable:/:Dynein','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTDynein','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttachActDynein')
react.VariableReferenceList = [['_', 'Variable:/:Dynein','-1']]
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTDynein','1']]
react.p = 0.001

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseKinesin')
diffuse.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
diffuse.D = 4e-12

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseDynein')
diffuse.VariableReferenceList = [['_', 'Variable:/:Dynein']]
diffuse.D = 4e-12

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:detach')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','1']]
react.SearchVacant = 1
react.k = 0.25

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:detachDynein')
react.VariableReferenceList = [['_', 'Variable:/:MTDyneinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Dynein','1']]
react.SearchVacant = 1
react.k = 0.25

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:hydrolysis')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.SearchVacant = 1
react.k = 100

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:hydrolysisDynein')
react.VariableReferenceList = [['_', 'Variable:/:MTDyneinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTDynein','1']]
react.SearchVacant = 1
react.k = 100

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:phosphorylate')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','1']]
react.SearchVacant = 1
react.k = 145

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:phosphorylateDynein')
react.VariableReferenceList = [['_', 'Variable:/:MTDynein','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTDyneinATP','1']]
react.SearchVacant = 1
react.k = 145

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:ratchet')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','1']]
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','1']]
react.BindingSite = 1
react.k = 55

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:ratchetDynein')
react.VariableReferenceList = [['_', 'Variable:/:MTDynein','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTDyneinATP','1']]
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','1']]
react.BindingSite = 0
react.k = 55

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:inactive')
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','1']]
react.k = 0

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePlus')
diffuse.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
diffuse.VariableReferenceList = [['_', 'Variable:/:actTubulin', '1']]
diffuse.D = 0.04e-12

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMTDynein')
diffuse.VariableReferenceList = [['_', 'Variable:/:MTDynein']]
diffuse.D = 0.04e-12

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:bindActDynein')
react.VariableReferenceList = [['_', 'Variable:/:MTDynein','-1']]
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','1']]
react.VariableReferenceList = [['_', 'Variable:/:MTDynein','1']]
react.ForcedSequence = 1
react.p = 1

#tagger = theSimulator.createEntity('TagProcess', 'Process:/:tagger')
#tagger.VariableReferenceList = [['_', 'Variable:/:GFP', '-1' ]]
#tagger.VariableReferenceList = [['_', 'Variable:/:Kinesin', '5' ]]
#tagger.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
#tagger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP']]

#life = theSimulator.createEntity('LifetimeLogProcess', 'Process:/:lifetime')
#life.VariableReferenceList = [['_', 'Variable:/:MTKinesin', '-1']]
#life.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP', '-1']]
#life.VariableReferenceList = [['_', 'Variable:/:Kinesin', '1']]

visualLogger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:visualLogger')
visualLogger.VariableReferenceList = [['_', 'Variable:/:Tubulin']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:TubulinM']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:TubulinP']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:Dynein']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:actTubulin']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTDynein' ]]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTDyneinATP' ]]
visualLogger.LogInterval = 1e-1

Microtubule = theSimulator.createEntity('MicrotubuleProcess', 'Process:/:Microtubule')
Microtubule.OriginX = 0
Microtubule.OriginY = 0
Microtubule.OriginZ = 0
Microtubule.RotateX = 0
Microtubule.RotateY = 0
Microtubule.RotateZ = RotateAngle
Microtubule.Radius = MTRadius
Microtubule.SubunitRadius = KinesinRadius
Microtubule.Length = MTLength
Microtubule.Filaments = 13
Microtubule.Periodic = 0
Microtubule.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
Microtubule.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
Microtubule.VariableReferenceList = [['_', 'Variable:/:MTDynein' ]]
Microtubule.VariableReferenceList = [['_', 'Variable:/:MTDyneinATP' ]]
Microtubule.VariableReferenceList = [['_', 'Variable:/:actTubulin' ]]
Microtubule.VariableReferenceList = [['_', 'Variable:/:Tubulin' , '-1']]
Microtubule.VariableReferenceList = [['_', 'Variable:/:TubulinM' , '-2']]
Microtubule.VariableReferenceList = [['_', 'Variable:/:TubulinP' , '-3']]

run(1000)

