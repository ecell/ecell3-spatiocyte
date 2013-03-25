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
theSimulator.createEntity('Variable', 'Variable:/Surface:PIP3').Value = 1500
theSimulator.createEntity('Variable', 'Variable:/Surface:PIP3s').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PIP3_PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PIP3s_PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:PTEN').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PTEN').Value = 120

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
logger.VariableReferenceList = [['_', 'Variable:/:PTEN']]
#logger.VariableReferenceList = [['_', 'Variable:/:Vacant']]
#logger.VariableReferenceList = [['_', 'Variable:/:Interface']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP3']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP2']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s']]
logger.LogInterval = 1e-5
logger.MultiscaleStructure = 0

#iterator = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iterate')
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN']]
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN']]
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:PIP3']]
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s']]
#iterator.LogInterval = 1e-5
#iterator.LogEnd = 0.009

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:PIP3']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PIP2']]
populator.VariableReferenceList = [['_', 'Variable:/:PTEN']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN']]
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

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactPTEN_PIP3')
react.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '1']]
react.p = 1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactPTEN_PIP2')
react.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '1']]
react.p = 1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:dissocPIP3_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Lipid', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '1']]
react.p = 0.3

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:dissocPIP2_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Lipid', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '1']]
react.p = 0.3

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:nucleatePIP3s')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:nucleatePIP2s')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:nucleatePIP23s')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:nucleatePIP3s_PIP3s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:nucleatePIP2s_PIP2s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:nucleatePIP3s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:nucleatePIP2s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:nucleatePIP23s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:extendPIP3s')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:extendPIP2s')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:extendPIP23s')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:extendPIP32s')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPIP3s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPIP2s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPIP3s_PIP3s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPIP2s_PIP2s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPIP2s_PIP3s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPIP3s_PIP2s_PTEN')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPIP3s_PIP3s_PTEN2')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPIP2s_PIP2s_PTEN2')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPIP2s_PIP3s_PTEN2')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPIP3s_PIP2s_PTEN2')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN', '1']]
react.p = 0.1

#react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:desorp')
#react.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/:PTEN', '1']]
#react.k = 20000

#react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissocPIP3s')
#react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '0']]
#react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '1']]
#react.Deoligomerize = 0
#react.ImplicitUnbind = 1
#react.SearchVacant = 0
#react.k = 1.5e+6

#react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissocPIP3sLip')
react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocPIP3sLip')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '1']]
react.Deoligomerize = 6
react.SearchVacant = 1
react.k = 1e+5

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocPIP2sLip')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP2', '1']]
react.Deoligomerize = 6
react.SearchVacant = 1
react.k = 1e+5

#react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissocPTENPIP3s')
react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocPTENPIP3s')
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '1']]
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

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP3')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP3']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP2')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP2']]
diffuser.D = 10e-12


diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP3_PTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP2_PTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP2_PTEN']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP3s')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP3s_PTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP2s')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP2s_PTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PIP2s_PTEN']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
diffuser.D = 0

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiPTEN_PIP3s')
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '1']]

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiPTEN_PIP3')
multi.VariableReferenceList = [['_', 'Variable:/Surface:PTEN', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '1']]
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
iterator.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
iterator.Iterations = 1
iterator.LogEnd = 0.0099
iterator.LogInterval = 1e-3
iterator.FrameDisplacement = 1

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/Surface:PTEN']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PIP3_PTEN', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s_PTEN', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PIP3', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PIP3s', '1']]
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


