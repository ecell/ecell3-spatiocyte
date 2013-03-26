sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 2.74e-9
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 1.0e-8
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 1e-7
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 1e-7
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 3

# Create the surface compartment:
theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')

theSimulator.createEntity('Variable', 'Variable:/Surface:PIP2').Value = 1500
theSimulator.createEntity('Variable', 'Variable:/Surface:PIP2s').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PIP2_PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PIP2s_PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:ANIO').Value = 1500 
theSimulator.createEntity('Variable', 'Variable:/Surface:ANIOs').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:ANIO_PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:ANIOs_PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PTEN').Value = 120 

#logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
#logger.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
#logger.VariableReferenceList = [['_', 'Variable:/:PTEN']]
#logger.VariableReferenceList = [['_', 'Variable:/:Vacant']]
#logger.VariableReferenceList = [['_', 'Variable:/:Interface']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:ANIO']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP2']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s']]
#logger.LogInterval = 1e-5
#logger.MultiscaleStructure = 0

#iterator = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iterate')
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN']]
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN']]
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:ANIO']]
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs']]
#iterator.LogInterval = 1e-5
#iterator.LogEnd = 0.009

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:ANIO']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PIP2']]
populator.VariableReferenceList = [['_', 'Variable:/:PTEN']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN']]

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

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactPTEN_PIP2')
react.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '1']]
react.p = 1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:dissocANIO_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Lipid', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '1']]
react.p = 0.3

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:dissocPIP2_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Lipid', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '1']]
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

#AC: ANIOx + PIP2x => 3 reactions
react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reactAC')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reactACs')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reactAsC')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

#AD: ANIOx + PIP2x_PTEN => 3 reactions
react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactAD')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactADs')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactAsD')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
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

#BC: ANIOx_PTEN + PIP2x => 3 reactions
react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactBC')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactBCs')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactBsC')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

#BD: ANIOx_PTEN + PIP2x_PTEN => 3 reactions
react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactBD')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactBDs')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactBsD')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

#CC: PIP2x + PIP2x => 2 reactions
react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reactCC')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reactCCs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

#CD: PIP2x + PIP2x_PTEN => 3 reactions
react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactCD')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactCDs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactCsD')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

#DD: PIP2x_PTEN + PIP2x_PTEN => 2 reactions
react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactDD')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactDDs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

#First order reactions (deoligomerizations)
react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocANIOsLip')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '1']]
react.Deoligomerize = 6
react.SearchVacant = 1
react.k = 1e+5

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocPIP2sLip')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '1']]
react.Deoligomerize = 6
react.SearchVacant = 1
react.k = 1e+5

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocPTENANIOs')
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '1']]
react.Deoligomerize = 6
react.SearchVacant = 1
react.k = 1e+5

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocPTENPIP2s')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '1']]
react.Deoligomerize = 6
react.SearchVacant = 1
react.k = 1e+5

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
diffuser.Interval = 6e-7
diffuser.Propensity = 1
diffuser.Origins = 1

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseANIO')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:ANIO']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP2')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP2']]
diffuser.D = 10e-12


diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseANIO_PTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP2_PTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseANIOs')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseANIOs_PTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP2s')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP2s_PTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN']]
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

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiPTEN_PIP2s')
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '1']]

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiPTEN_PIP2')
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '1']]

iterator = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iterate')
iterator.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s']]
iterator.VariableReferenceList = [['_', 'Variable:/Surface:PIP2']]
iterator.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs']]
iterator.VariableReferenceList = [['_', 'Variable:/Surface:ANIO']]
iterator.Iterations = 200
iterator.LogEnd = 0.0099
iterator.LogInterval = 1e-5
iterator.FileName = "anioPip2.csv"

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:ANIO_PTEN', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs_PTEN', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:ANIO', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:ANIOs', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
fil.Length = 1e-7
fil.Width = 1e-7
#fil.Filaments = 4
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
run(0.01)
end = time.time()
duration = end-start
print duration
