sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 2.74e-9
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 1.0e-8
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 1.3e-7
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 1.3e-7
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 3

# Create the surface compartment:
theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')

theSimulator.createEntity('Variable', 'Variable:/Surface:PG').Value = 3000
theSimulator.createEntity('Variable', 'Variable:/Surface:PGs').Value = 1
theSimulator.createEntity('Variable', 'Variable:/Surface:PG_MinD').Value = 128
theSimulator.createEntity('Variable', 'Variable:/Surface:PGs_MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:MinD').Value = 25

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
logger.VariableReferenceList = [['_', 'Variable:/:MinD']]
#logger.VariableReferenceList = [['_', 'Variable:/:Vacant']]
#logger.VariableReferenceList = [['_', 'Variable:/:Interface']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
logger.LogInterval = 0.00001

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
populator.VariableReferenceList = [['_', 'Variable:/:MinD']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop2')
populator.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
populator.UniformRadiusY = 0.99
populator.UniformRadiusZ = 0.99

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:adsorp')
react.VariableReferenceList = [['_', 'Variable:/:MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Vacant', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:react')
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:dissocPG_MinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Lipid', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
react.p = 0

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:dimerizePG')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 1

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:desorp')
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:MinD', '1']]
react.k = 20000

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissocPGs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
react.k = 10000

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMinDv')
diffuser.VariableReferenceList = [['_', 'Variable:/:MinD']]
diffuser.D = 8e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
diffuser.D = 3e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:propenMinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
diffuser.D = 4e-12
diffuser.Propensity = 1

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePG')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
diffuser.D = 5e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePGs')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePG_MinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
diffuser.D = 5e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePGs_MinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
diffuser.D = 0

multi = theSimulator.createEntity('MultiscaleProcess', 'Process:/:multiA')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]

multi = theSimulator.createEntity('MultiscaleProcess', 'Process:/:multiB')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]

multi = theSimulator.createEntity('MultiscaleProcess', 'Process:/:multiC')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]

multi = theSimulator.createEntity('MultiscaleProcess', 'Process:/:multiD')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
fil.Length = 1.3e-7
fil.Width = 1.3e-7
#fil.Filaments = 4
fil.SubunitRadius = 2.74e-9
fil.DiffuseRadius = 0.436e-9
fil.LipidRadius = 0.436e-9
fil.Periodic = 0

run(0.1)

