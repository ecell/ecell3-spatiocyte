sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 4.4e-9 
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 5e-8
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 1e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 1e-6
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 4

# Create the surface compartment:
theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')
theSimulator.createEntity('Variable', 'Variable:/Surface:MinD').Value = 300 
theSimulator.createEntity('Variable', 'Variable:/:MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PG').Value = 3112
theSimulator.createEntity('Variable', 'Variable:/Surface:PGs').Value = 500

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/:MinD'], ['_', 'Variable:/Surface:MinD'], ['_', 'Variable:/Surface:PG'], ['_', 'Variable:/Surface:PGs']]
logger.LogInterval = 0.0001

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]

binder = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:PGs2PG')
binder.VariableReferenceList = [['_', 'Variable:/Surface:PGs','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:PG','1']]
binder.k = 300


binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction0')
binder.VariableReferenceList = [['_', 'Variable:/:MinD','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:PG','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:MinD','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction2')
binder.VariableReferenceList = [['_', 'Variable:/:MinD','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:PGs','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:MinD','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction1')
binder.VariableReferenceList = [['_', 'Variable:/Surface:MinD','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:PG','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:PGs','1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:MinD','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:dissociate')
binder.VariableReferenceList = [['_', 'Variable:/Surface:MinD','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:VACANT','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:PGs','1']]
binder.VariableReferenceList = [['_', 'Variable:/:MinD','1']]
binder.p = 0.01

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
diffuser.D = 1e-13

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePG')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
diffuser.D = 1e-13

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePGs')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
diffuser.D = 1e-15

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMinDc')
diffuser.VariableReferenceList = [['_', 'Variable:/:MinD']]
diffuser.D = 1e-12

run(1)
