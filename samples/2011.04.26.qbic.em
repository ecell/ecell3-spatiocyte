Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 6e-8;    # m
}

System System(/)
{
  StepperID       SS; 
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
      Value 30e-6;      # m
    } 
  Variable Variable(LENGTHY)
    {
      Value 21e-6;      # m
    } 
  Variable Variable(LENGTHZ)
    {
      Value 21e-6;      # m
    } 
  Variable Variable(VACANT)
    {
      Value 0; 
    } 
  Process VisualizationLogProcess(loggerMean)
    {
      VariableReferenceList [_ Variable:/Sphere/Surface:A]
                            [_ Variable:/Sphere:A]
                            [_ Variable:/Sphere/Surface:LIPID]
                            [_ Variable:/Rod/Surface:LIPID];
      LogInterval 1;
    }
  Process MoleculePopulateProcess(populate)
    {
      VariableReferenceList [_ Variable:/Sphere/Surface:A]
                            [_ Variable:/Sphere:A];
    }
}

System System(/Sphere)
{
  StepperID       SS; 
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
      Value 20e-6;      # m
    } 
  Variable Variable(LENGTHY)
    {
      Value 20e-6;      # m
    } 
  Variable Variable(LENGTHZ)
    {
      Value 20e-6;      # m
    } 
  Variable Variable(VACANT)
    {
      Value 0; 
    } 
  Variable Variable(A)
    {
      Value 1000;      
    } 
  Process DiffusionProcess(diffuseA)
    {
      VariableReferenceList [_ Variable:.:A];
      D 0.2e-12;
    }
}


System System(/Sphere/Surface)
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
  Variable Variable(A)
    {
      Value 1000; 
    } 
  Process DiffusionProcess(diffuseA)
    {
      VariableReferenceList [_ Variable:.:A];
      D 0.2e-12;
    }
}



System System(/Rod)
{
  StepperID       SS; 
  Variable Variable(SHAPE)
    {
      Value 1;         # { 0: Spherical (uses SIZE) 
                       #   1: Rod (uses SIZE, LENGTHY=2*radius)
                       #   2: Cubic (uses SIZE)
                       #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                       #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable(LENGTHY)
    {
      Value 6e-6;      # m
    } 
  Variable Variable(ORIGINX)
    {
      Value 0.69;      # m
    } 
  Variable Variable(SIZE)
    {
      Value 160e-18;     # m^3
    } 
  Variable Variable(VACANT)
    {
      Value 1;         # { 0: full occupancy of volume
                       #   1: cut off intersecting volume }
    } 
  Variable Variable(DIFFUSIVE)
    {
      Name "/:Sphere";
    }
}


System System(/Rod/Surface)
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
  Variable Variable(DIFFUSIVE)
    {
      Name "/Sphere:Surface";
    }
}
