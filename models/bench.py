sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 2.74e-9
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 1e-8
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

theSimulator.createEntity('Variable', 'Variable:/Surface:PG').Value = 8000
theSimulator.createEntity('Variable', 'Variable:/Surface:PGs').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PG_MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PGs_MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:MinD').Value = 200

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
logger.VariableReferenceList = [['_', 'Variable:/:MinD']]
#logger.VariableReferenceList = [['_', 'Variable:/:Vacant']]
#logger.VariableReferenceList = [['_', 'Variable:/:Lipid']]
#logger.VariableReferenceList = [['_', 'Variable:/:Interface']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
#logger.VariableReferenceList = [['_', 'Variable:/Surface:VACANT']]
logger.LogInterval = 1e-3

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
populator.VariableReferenceList = [['_', 'Variable:/:MinD']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop22')
populator.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
#populator.UniformRadiusY = 0.1
#populator.UniformRadiusZ = 0.1

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop2')
populator.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]

#diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMinDv')
#diffuser.VariableReferenceList = [['_', 'Variable:/:MinD']]
#diffuser.D = 5e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
diffuser.D = 5e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePG')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePG_MinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
diffuser.D = 0

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiA')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiB')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiC')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiD')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
#fil.Length = 1e-7
#fil.Width = 1e-7
#fil.Autofit = 1
#fil.Filaments = 1
#fil.Subunits = 1
fil.SubunitRadius = 1.74e-9
fil.DiffuseRadius = 0.436e-9
fil.LipidRadius = 0.436e-9
fil.RegularLattice = 1
#fil.DiffuseRadius = 1.74e-9
#fil.LipidRadius = 1.74e-9
fil.Periodic = 1

#run(0.0001)
import time
run(1e-6)
print "Done stirring. Now running..."
start = time.time()
run(0.0001)
end = time.time()
duration = end-start
print duration

