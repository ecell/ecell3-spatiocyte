# Reaction between molecules from distinct surfaces
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 6e-8;    # m
}

System System(/)
{
  StepperID       SS; 
  Variable Variable(TYPE)
    {
      Value 0;         # { 0: Volume
                       #   1: Surface }
    } 
  Variable Variable(SHAPE)
    {
      Value 3;         # { 0: Spherical (uses SIZE) 
                       #   1: Rod (uses SIZE, LENGTHY=2*radius)
                       #   2: Cubic (uses SIZE)
                       #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                       #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable(LENGTHX)
    {
      Value 8e-6;      # m
    } 
  Variable Variable(LENGTHY)
    {
      Value 8e-6;      # m
    } 
  Variable Variable(LENGTHZ)
    {
      Value 1e-6;      # m
    } 
  Variable Variable(XYPLANE)
    {
      Value 3;         # { 0: REFLECTIVE
                       #   1: PERIODIC
                       #   2: UNIPERIODIC
                       #   3: REMOVE_UPPER
                       #   4: REMOVE_LOWER
                       #   5: REMOVE_BOTH }
    }
  Variable Variable(XZPLANE)
    {
      Value 5;         # { 0: REFLECTIVE
                       #   1: PERIODIC
                       #   2: UNIPERIODIC
                       #   3: REMOVE_UPPER
                       #   4: REMOVE_LOWER
                       #   5: REMOVE_BOTH }
    } 
  Variable Variable(YZPLANE)
    {
      Value 5;         # { 0: REFLECTIVE
                       #   1: PERIODIC
                       #   2: UNIPERIODIC
                       #   3: REMOVE_UPPER
                       #   4: REMOVE_LOWER
                       #   5: REMOVE_BOTH }
    } 
  Variable Variable(VACANT)
    {
      Value 0; 
    } 
  Process VisualizationLogProcess(loggerMean)
    {
      VariableReferenceList [_ Variable:/Membrane:A]
                            [_ Variable:/Cell/Membrane:B]
                            [_ Variable:/Cell/Membrane:C]
                            [_ Variable:/Membrane:LIPID]
                            [_ Variable:/Cell/Membrane:LIPID];
    }
  Process MoleculePopulateProcess(populateUniform)
    {
      VariableReferenceList [_ Variable:/Membrane:A]
                            [_ Variable:/Cell/Membrane:B]
                            [_ Variable:/Cell/Membrane:C];
    }
}

System System(/Membrane)
{
  StepperID SS;

  Variable Variable(TYPE)
    {
      Value 1;         # { 0: Volume
                       #   1: Surface }
    } 
  Variable Variable(LIPID)
    {
      Value 0;         # { 0: open surface (cut off overlapped parts)
                       #   1: enclosed surface }
    } 
  Variable Variable(REACTIVE)
    {
      Name "/Cell:Membrane";
    }
  Variable Variable(A)
    {
      Value 100;       # molecule number 
    } 
  Process DiffusionProcess(diffuseMembranePTEN)
    {
      VariableReferenceList [_ Variable:/Membrane:A];
      D 0.2e-12;       # m^2/s
    }
}



System System(/Cell)
{
  StepperID       SS; 
  Variable Variable(TYPE)
    {
      Value 0;         # { 0: Volume
                       #   1: Surface }
    } 
  Variable Variable(SHAPE)
    {
      Value 4;         # { 0: Spherical (uses SIZE) 
                       #   1: Rod (uses SIZE, LENGTHY=2*radius)
                       #   2: Cubic (uses SIZE)
                       #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                       #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable(LENGTHX)
    {
      Value 6e-6;      # m
    } 
  Variable Variable(LENGTHY)
    {
      Value 6e-6;      # m
    } 
  Variable Variable(LENGTHZ)
    {
      Value 0.8e-6;    # m
    } 
  Variable Variable(ORIGINZ)
    {
      Value -0.5;       
    } 
  Variable Variable(VACANT)
    {
      Value 0; 
    } 
}

System System(/Cell/Membrane)
{
  StepperID SS;

  Variable Variable(TYPE)
    {
      Value 1;         # { 0: Volume
                       #   1: Surface }
    } 
  Variable Variable(LIPID)
    {
      Value 1;         # { 0: open surface (cut off overlapped parts)
                       #   1: enclosed surface }
    } 
  Variable Variable(B)
    {
      Value 50;        # molecule number 
    } 
  Variable Variable(C)
    {
      Value 0;         # molecule number
    } 
  Process DiffusionProcess(diffusePI3K_PIP2)
    {
      VariableReferenceList [_ Variable:/Cell/Membrane:B];
      D 0.2e-12;       # m^2/s
    }

  Process DiffusionProcess(diffusePTEN_PIP2)
    {
      VariableReferenceList [_ Variable:/Cell/Membrane:C];
      D 0.2e-12;       # m^2/s
    }
  Process DiffusionInfluencedReactionProcess(phosphataseBind)
    {
      VariableReferenceList [_ Variable:/Membrane:A -1]
                            [_ Variable:/Cell/Membrane:B -1]
                            [_ Variable:/Cell/Membrane:C 1];
      p 1;
    }
}

