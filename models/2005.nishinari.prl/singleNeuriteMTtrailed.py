import math
import random

Protofilaments = 13
RotateAngle = 0 #math.pi/4
MTRadius = 12.5e-9
VoxelRadius = 0.4e-8
KinesinRadius = 0.4e-8
minDist = 75e-9
dendriteRadius = 0.75e-6
dendriteLength = 10e-6
totalMTLength = 620e-6/2

lengths = [0.1e-6, 4e-6, 6e-6, 8e-6, 9.91e-6]
lengthFreqs = [121, 13, 12, 9]
maxMTLength = lengths[len(lengths)-1]

maxBoundaryMT = 0.95

maxLengthCnt = 0
while(totalMTLength-maxMTLength > 0):
  maxLengthCnt = maxLengthCnt+1
  totalMTLength = totalMTLength-maxMTLength
remainderLength = totalMTLength

def getLengthIndex(length):
  for i in range(len(lengths)):
    if(length < lengths[i]):
      return i-1

def getShortestLength():
  for i in range(len(lengthFreqs)):
    if(lengthFreqs[i]):
      return lengths[i+1]

mtOriginX = []
mtOriginZ = []
mtOriginY = []
completedLengths = []

def isNotInside(currX, sign, currLength):
  x = currX-sign*(currLength/2)/(dendriteLength/2)
  maxX = x*dendriteLength/2 + currLength/2
  if(maxX > dendriteLength*maxBoundaryMT/2):
    return True
  minX = x*dendriteLength/2 - currLength/2
  if(minX < -dendriteLength*maxBoundaryMT/2):
    return True
  return False

def isSpacedOut(x, y, z, length):
  for i in range(len(completedLengths)):
    maxOriX = mtOriginX[i]*dendriteLength/2 + completedLengths[i]/2
    minOriX = mtOriginX[i]*dendriteLength/2 - completedLengths[i]/2
    maxX = x*dendriteLength/2 + length/2
    minX = x*dendriteLength/2 - length/2
    y2 = math.pow((y-mtOriginY[i])*dendriteRadius, 2)
    z2 = math.pow((z-mtOriginZ[i])*dendriteRadius, 2)
    if((minX <= maxOriX or maxX >= minOriX) and math.sqrt(y2+z2) < minDist):
      return False
    elif(minX > maxOriX and math.sqrt(y2+z2+math.pow(minX-maxOriX, 2)) < minDist):
      return False
    elif(maxX < minOriX and math.sqrt(y2+z2+math.pow(maxX-minOriX, 2)) < minDist):
      return False
  return True

def updateOrigin(currX, sign, currLength):
  x = currX-sign*(currLength/2)/(dendriteLength/2)
  y = random.uniform(-maxBoundaryMT, maxBoundaryMT)
  z = random.uniform(-maxBoundaryMT, maxBoundaryMT)
  while(y*y+z*z > 0.9 or not isSpacedOut(x, y, z, currLength)):
    y = random.uniform(-maxBoundaryMT, maxBoundaryMT)
    z = random.uniform(-maxBoundaryMT, maxBoundaryMT)
  mtOriginX.append(x)
  mtOriginY.append(y)
  mtOriginZ.append(z)
  completedLengths.append(currLength)

sign = 1
x = maxBoundaryMT
for i in xrange(maxLengthCnt):
  if(sign == 1):
    sign = -1
  else:
    sign = 1
  x = maxBoundaryMT*sign
  availableLength = maxMTLength 
  while(availableLength > getShortestLength()):
    length = random.uniform(lengths[0], lengths[len(lengths)-1])
    index = getLengthIndex(length)
    while(length > availableLength or lengthFreqs[index] == 0 or isNotInside(x, sign, length)):
      length = random.uniform(lengths[0], lengths[len(lengths)-1])
      index = getLengthIndex(length)
    availableLength = availableLength-length
    lengthFreqs[index] = lengthFreqs[index]-1
    updateOrigin(x, sign, length)
    x = x-sign*length/(dendriteLength/2)

shuffledIndices = []
for i in range(len(lengthFreqs)):
  for j in range(int(lengthFreqs[i])):
    shuffledIndices.append(i)
random.shuffle(shuffledIndices)

sign = 1
x = maxBoundaryMT
while(len(shuffledIndices)):
  if(sign == 1):
    sign = -1
  else:
    sign = 1
  x = maxBoundaryMT*sign
  tmp = []
  for i in range(len(shuffledIndices)):
    index = shuffledIndices[i]
    currLength = random.uniform(lengths[index], lengths[index+1])
    if(not isNotInside(x, sign, currLength)):
      oriX = x
      x = x-sign*(currLength/2)/(dendriteLength/2)
      y = random.uniform(-maxBoundaryMT, maxBoundaryMT)
      z = random.uniform(-maxBoundaryMT, maxBoundaryMT)
      while(y*y+z*z > 0.9 or not isSpacedOut(x, y, z, lengths[index])):
        y = random.uniform(-maxBoundaryMT, maxBoundaryMT)
        z = random.uniform(-maxBoundaryMT, maxBoundaryMT)
      mtOriginX.append(x)
      mtOriginY.append(y)
      mtOriginZ.append(z)
      completedLengths.append(currLength)
      x = x-sign*(currLength/2)/(dendriteLength/2)
    else:
      tmp.append(index)
  shuffledIndices = tmp

theSimulator.createStepper('SpatiocyteStepper', 'SS').VoxelRadius = VoxelRadius
theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 2
theSimulator.createEntity('Variable', 'Variable:/:ROTATEZ').Value = RotateAngle
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = dendriteLength
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = dendriteRadius*2
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:Kinesin').Value = 1000
theSimulator.createEntity('Variable', 'Variable:/:MTKinesin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MTKinesinATP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Tubulin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:actTubulin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:TubulinM' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:TubulinP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:GFP' ).Value = 0

theSimulator.createEntity('System', 'System:/:Membrane').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Membrane:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Membrane:VACANT')

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populate')
populator.VariableReferenceList = [['_', 'Variable:/:Kinesin']]

#tagger = theSimulator.createEntity('TagProcess', 'Process:/:tagger')
#tagger.VariableReferenceList = [['_', 'Variable:/:GFP', '-1' ]]
#tagger.VariableReferenceList = [['_', 'Variable:/:Kinesin', '20' ]]
#tagger.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
#tagger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP']]

visualLogger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:visualLogger')
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesin', '10600']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP', '10600' ]]
visualLogger.VariableReferenceList = [['_', 'Variable:/:actTubulin', '10400' ]]
visualLogger.VariableReferenceList = [['_', 'Variable:/:Tubulin' ]]
visualLogger.LogInterval = 0.5
#visualLogger.FileName = "tmp.dat"

#life = theSimulator.createEntity('LifetimeLogProcess', 'Process:/:lifetime')
#life.VariableReferenceList = [['_', 'Variable:/:MTKinesin', '-1']]
#life.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP', '-1']]
#life.VariableReferenceList = [['_', 'Variable:/:Kinesin', '1']]
#life.Iterations = 500
#life.LogEnd = 99
#life.FileName = "LifetimeLogKon.csv"

#iterate = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iter')
#iterate.VariableReferenceList = [['_', 'Variable:/:Kinesin', '-1']]
#iterate.VariableReferenceList = [['_', 'Variable:/:Tubulin', '-1']]
#iterate.Iterations = 1
#iterate.LogInterval = 1e-2
#iterate.SaveCounts = 100
#iterate.FileName = "IterateLog293792p0.1268.csv"
#iterate.LogEnd = 30
#iterate.LogStart = 1

#micro = theSimulator.createEntity('MicroscopyTrackingProcess', 'Process:/:track')
#micro.VariableReferenceList = [['_', 'Variable:/:Tubulin', '0']]
#micro.VariableReferenceList = [['_', 'Variable:/:actTubulin', '1']]
#micro.VariableReferenceList = [['_', 'Variable:/:MTKinesin', '1']]
#micro.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP', '1']]
#micro.VariableReferenceList = [['_', 'Variable:/:actTubulin', '-1']]
#micro.VariableReferenceList = [['_', 'Variable:/:MTKinesin', '-1']]
#micro.VariableReferenceList = [['_', 'Variable:/:MTKinesin', '-1']]
#micro.VariableReferenceList = [['_', 'Variable:/:MTKinesin', '12000']]
#micro.VariableReferenceList = [['_', 'Variable:/:actTubulin', '10700']]
#micro.VariableReferenceList = [['_', 'Variable:/:Tubulin', '10300']]
#micro.ExposureTime = 0.1
#micro.LogInterval = 0.01
#micro.FileName = "microKinesin_actD4_kon2.586e-22-23_koff0.25_m100_r20.dat"

#clogger = theSimulator.createEntity('CoordinateLogProcess', 'Process:/:clog')
#clogger.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
#clogger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
#clogger.LogInterval = 0.03
#clogger.FileName = "CoordinateLog_act0.15_m100.csv"

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:detachPlus')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinP','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttach')
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.k = 2.5863133e-23

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttachAct')
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.k = 2.5863133e-21

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseKinesin')
diffuse.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
diffuse.D = 4e-12

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:detach')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','1']]
react.SearchVacant = 1
react.k = 0.25

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:hydrolysis')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.SearchVacant = 1
react.k = 100

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:phosphorylate')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','1']]
react.SearchVacant = 1
react.k = 145

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:ratchet')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','1']]
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','1']]
react.BindingSite = 1
react.k = 55

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:inactive')
react.VariableReferenceList = [['_', 'Variable:/:actTubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','1']]
react.k = 2

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePlus')
diffuse.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
diffuse.VariableReferenceList = [['_', 'Variable:/:actTubulin', '1']]
diffuse.D = 0.04e-12

#Histogram = theSimulator.createEntity('HistogramLogProcess', 'Process:/:Histogram')
#Histogram.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
#Histogram.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
#Histogram.VariableReferenceList = [['_', 'Variable:/:Kinesin' ]]
#Histogram.RotateZ = RotateAngle
#Histogram.Length = dendriteLength
#Histogram.Radius = dendriteRadius
#Histogram.Bins = 20 
#Histogram.LogInterval = 1
#Histogram.FileName = "SingleNeurite.HisLog1.00_25s.csv"
#Histogram.LogEnd = 25
#Histogram.Iterations = 50

for i in range(len(completedLengths)):
  Microtubule = theSimulator.createEntity('MicrotubuleProcess', 'Process:/:Microtubule%d' %i)
  Microtubule.OriginX = mtOriginX[i]
  Microtubule.OriginY = mtOriginY[i]
  Microtubule.OriginZ = mtOriginZ[i]
  Microtubule.RotateX = 0
  Microtubule.RotateY = 0
  Microtubule.RotateZ = RotateAngle
  Microtubule.Length = completedLengths[i]
  Microtubule.Radius = MTRadius
  Microtubule.SubunitRadius = KinesinRadius
  Microtubule.Filaments = Protofilaments
  Microtubule.Periodic = 0
  Microtubule.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
  Microtubule.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
  Microtubule.VariableReferenceList = [['_', 'Variable:/:actTubulin']]
  Microtubule.VariableReferenceList = [['_', 'Variable:/:Tubulin' , '-1']]
  Microtubule.VariableReferenceList = [['_', 'Variable:/:TubulinM' , '-2']]
  Microtubule.VariableReferenceList = [['_', 'Variable:/:TubulinP' , '-3']]
run(300.1)




