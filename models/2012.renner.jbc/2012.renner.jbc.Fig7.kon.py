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

theSimulator.createEntity('Variable', 'Variable:/Surface:PG').Value = 15048
theSimulator.createEntity('Variable', 'Variable:/Surface:PGs').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PG_MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PGs_MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:MinD').Value = 0

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
logger.LogInterval = 1e-5
logger.MultiscaleStructure = 1

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
populator.VariableReferenceList = [['_', 'Variable:/:MinD']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop2')
populator.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
populator.UniformRadiusY = 0.99
populator.UniformRadiusZ = 0.99

#react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:adsorp')
#react.VariableReferenceList = [['_', 'Variable:/:MinD', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/:Vacant', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '1']]
#react.p = 1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactMinD_PG')
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
react.p = 1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:dissocPG_MinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Lipid', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
react.p = 0.3


#AA: PGx + PGx => 2 reactions
react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reactAA')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reactAAs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.p = 0.1

#AB: PGx + PGx_MinD => 3 reactions
react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactAB')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactABs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactAsB')
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 0.1


#BB: PGx_MinD + PGx_MinD => 2 reactions
react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactBB')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:reactBBs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 0.1

#First order reactions (deoligomerizations)
react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocPGsLip')
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
react.Deoligomerize = 6
react.SearchVacant = 1
react.k = 5.5e+5

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:dissocMinDPGs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
react.Deoligomerize = 6
react.SearchVacant = 1
react.k = 5.5e+5

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:konMinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '1']]
react.k = 2.92e+20

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:koffMinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Vacant', '1']]
react.k = 1e+6

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:propenMinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
diffuser.Interval = 9e-8
diffuser.Propensity = 1.7
#diffuser.Propensity = 1
#diffuser.Origins = 1

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePG')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePG_MinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePGs')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePGs_MinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
diffuser.D = 0

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiMinD_PGs')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '1']]

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiMinD_PG')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '1']]

#iterator = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iterate')
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]
#iterator.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
#iterator.Iterations = 100
#iterator.LogEnd = 0.0099
#iterator.LogInterval = 1e-5
#iterator.SaveCounts = 20
#iterator.FileName = "stpl.csv"

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
fil.SubunitRadius = 2.304e-9
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
