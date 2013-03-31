sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 2.74e-9
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 1.0e-8
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 0.3e-7
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 0.3e-7
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 3

# Create the surface compartment:
theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')

theSimulator.createEntity('Variable', 'Variable:/Surface:ANIO').Value = 146
theSimulator.createEntity('Variable', 'Variable:/Surface:ANIOs').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:ANIO_PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:ANIOs_PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PTEN').Value = 1

#logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
#logger.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:ANIO']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs']]
#logger.LogInterval = 1e+5
#logger.MultiscaleStructure = 0

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:ANIO']]
populator.VariableReferenceList = [['_', 'Variable:/:PTEN']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop2')
populator.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
populator.UniformRadiusY = 0.99
populator.UniformRadiusZ = 0.99

#react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:adsorp')
#react.VariableReferenceList = [['_', 'Variable:/:PTEN', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/:Vacant', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '1']]
#react.p = 1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactPTEN_ANIO')
react.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '1']]
react.p = 1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:dissocANIO_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Lipid', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '1']]
react.p = 0.3


#AA: ANIOx + ANIOx => 2 reactions
react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reactAA')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reactAAs')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.p = 0.1

#AB: ANIOx + ANIOx_PTEN => 3 reactions
react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactAB')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactABs')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactAsB')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.p = 0.1


#BB: ANIOx_PTEN + ANIOx_PTEN => 2 reactions
react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactBB')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactBBs')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.p = 0.1

#First order reactions (deoligomerizations)
react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocANIOsLip')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '1']]
react.Deoligomerize = 6
react.SearchVacant = 1
react.k = 5e+5

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocPTENANIOs')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '1']]
react.Deoligomerize = 6
react.SearchVacant = 1
react.k = 5e+5

#diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePTENv')
#diffuser.VariableReferenceList = [['_', 'Variable:/:PTEN']]
#diffuser.D = 5e-12

#diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePTEN')
#diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
#diffuser.D = 5e-12

#rotator = theSimulator.createEntity('DiffusionProcess', 'Process:/:rotatePTEN')
#rotator.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
#rotator.D = 1e-12

#rotator = theSimulator.createEntity('DiffusionProcess', 'Process:/:rotatePropenPTEN')
#rotator.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
#rotator.D = 1e-12
#rotator.Propensity = 1

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:propenPTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
diffuser.Interval = 5e-7
diffuser.Propensity = 3.2
#diffuser.Propensity = 1
diffuser.Origins = 1

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseANIO')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:ANIO']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseANIO_PTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseANIOs')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseANIOs_PTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
diffuser.D = 0

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiPTEN_ANIOs')
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '1']]

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiPTEN_ANIO')
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '1']]


iterator = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iterate')
iterator.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
iterator.Iterations = 10000
iterator.SeparateFiles = 1
iterator.LogEnd = 5
iterator.LogStart = 0.006
iterator.LogInterval = 1e-3
iterator.FrameDisplacement = 1
iterator.FileStartCount = 300
iterator.FileName = "2013.matsuoka.pcb.Fig4B.sim.anionic.d5e5.p3.2csv"

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
fil.SubunitRadius = 1.74e-9
fil.SubunitAngle = 0
fil.DiffuseRadius = 0.436e-9
fil.LipidRadius = 0.436e-9
fil.Periodic = 1
fil.RegularLattice = 1

import time
run(1e-6)
print "Done stirring. Now running..."
start = time.time()
run(5.1)
end = time.time()
duration = end-start
print duration
