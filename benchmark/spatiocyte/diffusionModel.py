import time

sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 6e-8;

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 60e-6;
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 60e-6;
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 60e-6;
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:A').Value = 10000

#log = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:log')
#log.VariableReferenceList = [['_', 'Variable:.:A']]

pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
pop.VariableReferenceList = [['_', 'Variable:.:A']]

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diff')
dif.VariableReferenceList = [['_', 'Variable:.:A']]
dif.D = 1e-14

run(10)
start = time.time()
run(10000)
end = time.time()
print end-start
