sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 2.74e-9
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 1.0e-8
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 1.115e-7
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 1.115e-7
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 3

# Create the surface compartment:
theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')
theSimulator.createEntity('Variable', 'Variable:/Surface:MinD').Value = 0

#logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
#logger.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
#logger.LogInterval = 1
#logger.MultiscaleStructure = 1

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop2')
populator.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
populator.UniformRadiusY = 0.99
populator.UniformRadiusZ = 0.99

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:konMinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '1']]
react.k = 2.04e+15

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:koffMinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Vacant', '1']]
react.k = 0.0589

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:propenMinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
diffuser.D = 1e-15

iterator = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iterate')
iterator.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
iterator.Iterations = 100
iterator.LogEnd = 1000
iterator.LogStart = 100
iterator.LogInterval = 1

#life = theSimulator.createEntity('LifetimeLogProcess', 'Process:/:lifetime')
#life.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
#life.VariableReferenceList = [['_', 'Variable:/:Vacant', '1']]
#life.Iterations = 1
#life.LogEnd = 9.99
#life.FileName = "LifetimeLogMinDKon.csv"

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
fil.SubunitRadius = 2.81e-9
fil.SubunitAngle = 0
#fil.DiffuseRadius = 0.436e-9
fil.DiffuseRadius = 2.81e-9
#fil.LipidRadius = 0.436e-9
fil.LipidRadius = 2.81e-9
fil.Periodic = 1
fil.RegularLattice = 1

import time
run(1e-6)
print "Done stirring. Now running..."
start = time.time()
run(1001)
end = time.time()
duration = end-start
print duration
