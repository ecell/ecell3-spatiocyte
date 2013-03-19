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

theSimulator.createEntity('Variable', 'Variable:/Surface:PG').Value = 3000
theSimulator.createEntity('Variable', 'Variable:/Surface:PGs').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PG_MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:PGs_MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MinD').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:MinD').Value = 120

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
logger.VariableReferenceList = [['_', 'Variable:/:MinD']]
#logger.VariableReferenceList = [['_', 'Variable:/:Vacant']]
#logger.VariableReferenceList = [['_', 'Variable:/:Interface']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
logger.LogInterval = 1e-5
logger.MultiscaleStructure = 0

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
populator.VariableReferenceList = [['_', 'Variable:/:MinD']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop22')
populator.VariableReferenceList = [['_', 'Variable:/Surface:PGs']]
populator.UniformRadiusY = 0.1
populator.UniformRadiusZ = 0.1

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop2')
populator.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
populator.UniformRadiusY = 0.99
populator.UniformRadiusZ = 0.99

#react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:adsorp')
#react.VariableReferenceList = [['_', 'Variable:/:MinD', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/:Vacant', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '1']]
#react.p = 1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:react')
react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
react.p = 1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:dissocPG_MinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:Lipid', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
react.p = 0.4

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:nucleatePGs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:nucleatePGs_PGs_MinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 0.1

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:nucleatePGs_MinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 0.1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:extendPGs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.p = 0.5

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPGs_MinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 0.5

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPGs_PGs_MinD')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 0.5

react = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:extendPGs_PGs_MinD2')
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
react.p = 0.5

#react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:desorp')
#react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/:MinD', '1']]
#react.k = 20000

#react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissocPGs')
#react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
#react.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '0']]
#react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
#react.Deoligomerize = 0
#react.ImplicitUnbind = 1
#react.SearchVacant = 0
#react.k = 1.5e+6

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissocPGsLip')
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
react.Deoligomerize = 6
react.k = 216666

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissocMinDPGs')
react.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
react.Deoligomerize = 6
react.k = 216666

#diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMinDv')
#diffuser.VariableReferenceList = [['_', 'Variable:/:MinD']]
#diffuser.D = 5e-12

#diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseMinD')
#diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
#diffuser.D = 5e-12

rotator = theSimulator.createEntity('DiffusionProcess', 'Process:/:ro----------tateMinD')
rotator.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
rotator.D = 1e-12

#rotator = theSimulator.createEntity('DiffusionProcess', 'Process:/:rotatePropenMinD')
#rotator.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
#rotator.D = 1e-12
#rotator.Propensity = 1

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:pro--------------penMinD')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
diffuser.D = 10e-12
diffuser.Propensity = 1

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diff------------usePG')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:PG']]
diffuser.D = 10e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diff----------------------usePG_MinD')
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

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiC')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '1']]

multi = theSimulator.createEntity('MultiscaleReactionProcess', 'Process:/:multiD')
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG', '-1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
multi.VariableReferenceList = [['_', 'Variable:/Surface:MinD', '1']]

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/Surface:MinD']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PG_MinD', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PGs_MinD', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PG', '1']]
fil.VariableReferenceList = [['_', 'Variable:/Surface:PGs', '1']]
fil.Length = 1e-7
fil.Width = 1e-7
#fil.Filaments = 4
fil.SubunitRadius = 1.74e-9
fil.SubunitAngle = 0.7
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


