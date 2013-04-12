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
theSimulator.createEntity('Variable', 'Variable:/:Kinesin').Value = 8000
theSimulator.createEntity('Variable', 'Variable:/:MTKinesin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MTKinesinATP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Tubulin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:TubulinM' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:TubulinP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:GFP' ).Value = 0

theSimulator.createEntity('System', 'System:/:Membrane').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Membrane:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Membrane:VACANT')

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populate')
populator.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
populator.VariableReferenceList = [['_', 'Variable:/:GFP']]

#tagger = theSimulator.createEntity('TagProcess', 'Process:/:tagger')
#tagger.VariableReferenceList = [['_', 'Variable:/:GFP', '-1' ]]
#tagger.VariableReferenceList = [['_', 'Variable:/:Kinesin', '20' ]]
#tagger.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
#tagger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP']]

visualLogger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:visualLogger')
visualLogger.VariableReferenceList = [['_', 'Variable:/:GFP' ]]
visualLogger.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:Interface']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
visualLogger.LogInterval = 0.01

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttach')
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','0']]
react.p = 1

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseKinesin')
diffuse.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
diffuse.D = 1e-12

#microLogger = theSimulator.createEntity('MicroscopyTrackingProcess', 'Process:/:track')
#microLogger.VariableReferenceList = [['_', 'Variable:/:GFP', '1']]
#microLogger.VariableReferenceList = [['_', 'Variable:/:GFP', '-1']]


#Histogram = theSimulator.createEntity('HistogramLogProcess', 'Process:/:Histogram')
#Histogram.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
#Histogram.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
#Histogram.VariableReferenceList = [['_', 'Variable:/:Kinesin' ]]
#Histogram.RotateZ = RotateAngle
#Histogram.Length = dendriteLength
#Histogram.Radius = dendriteRadius
#Histogram.Bins = 20 
#Histogram.LogInterval = 1
#Histogram.FileName = "2012.08.17.singleNeurite.HisLog.csv"

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
  Microtubule.VariableReferenceList = [['_', 'Variable:/:Tubulin' , '-1']]
  Microtubule.VariableReferenceList = [['_', 'Variable:/:TubulinM' , '-2']]
  Microtubule.VariableReferenceList = [['_', 'Variable:/:TubulinP' , '-3']]
run(1000)




