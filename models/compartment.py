sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 2.74e-9
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 1.0e-8
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 1.0e-7
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 1.0e-7
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 3

# Create the surface compartment:
theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')

theSimulator.createEntity('Variable', 'Variable:/Surface:PG').Value = 3062
theSimulator.createEntity('Variable', 'Variable:/Surface:PG_MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:MinD').Value = 20

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
logger.VariableReferenceList = [['_', 'Variable:/:MinD']]
logger.VariableReferenceList = [['_', 'Variable:/:Vacant']]
#logger.VariableReferenceList = [['_', 'Variable:/:Interface']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
logger.LogInterval = 0.000001

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
populator.VariableReferenceList = [['_', 'Variable:/:MinD']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop2')
populator.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
populator.UniformRadiusY = 0.99
populator.UniformRadiusZ = 0.99

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:react')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Lipid', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
react.p = 0.1

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:desorp')
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:MinD', '1']]
react.k = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMinDv')
diffuser.VariableReferenceList = [['_', 'Variable:/:MinD']]
diffuser.D = 5e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
diffuser.D = 0.5e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePG')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
diffuser.D = 5e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePG_MinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
diffuser.D = 5e-12

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
fil.Length = 1.0e-7
fil.Width = 1.0e-7
#fil.Filaments = 4
fil.SubunitRadius = 2.74e-9
fil.DiffuseRadius = 0.436e-9
fil.LipidRadius = 0.436e-9
fil.Periodic = 1

run(0.0005)

