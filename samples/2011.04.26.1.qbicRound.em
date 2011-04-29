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
  Variable Variable(PIP2m)
    {
      Value             1970; 
    } 
  Variable Variable(PIP3m)
    {
      Value             0; 
    } 
  Variable Variable(PIP3a)
    {
      Value             0; 
    } 
  Variable Variable(PTENm)
    {
      Value             494; 
    } 
  Variable Variable(PI3Km)
    {
      Value             4940; 
    } 
  Variable Variable(PIP2)
    {
      Value             1.284e+4; 
      Name "HD";
    } 
  Variable Variable(PI3K)
    {
      Value             1.477e+4; 
      Name "HD";
    } 
  Variable Variable(PTEN)
    {
      Value             9880; 
      Name "HD";
    } 
  Process DiffusionProcess(diffusePIP2)
    {
      VariableReferenceList [_ Variable:.:PIP2m];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePIP3)
    {
      VariableReferenceList [_ Variable:.:PIP3m];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePIP3a)
    {
      VariableReferenceList [_ Variable:.:PIP3a];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePTEN)
    {
      VariableReferenceList [_ Variable:.:PTENm];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePI3K)
    {
      VariableReferenceList [_ Variable:.:PI3Km];
      D 1e-14;
    }

#Membrane recruitments:
  Process SpatiocyteNextReactionProcess(recruitPIP2)
    {
      VariableReferenceList   [_ Variable:.:PIP2 -1]
                              [_ Variable:.:PIP2m 1];
      k 4e-2;
    }
  Process SpatiocyteNextReactionProcess(recruitPTEN)
    {
      VariableReferenceList   [_ Variable:.:PTEN -1]
                              [_ Variable:.:PIP2m -1]
                              [_ Variable:.:PTENm 1]
                              [_ Variable:.:PIP2m 1];
      k 2e-14;
    }
  Process SpatiocyteNextReactionProcess(recruitPI3Ka)
    {
      VariableReferenceList   [_ Variable:.:PIP3a -1]
                              [_ Variable:.:PI3K -1]
                              [_ Variable:.:PIP3m 1]
                              [_ Variable:.:PI3Km 1];
      k 1e-13;
    }
#Membrane recruitments end


#Activations:
  Process DiffusionInfluencedReactionProcess(dimerPIP3)
    {
      VariableReferenceList   [_ Variable:.:PIP3m -1]
                              [_ Variable:.:PIP3m -1]
                              [_ Variable:.:PIP3a 1]
                              [_ Variable:.:PIP3a 1];
      p 0.65;
    }
  Process DiffusionInfluencedReactionProcess(PIP2toPIP3)
    {
      VariableReferenceList   [_ Variable:.:PIP2m -1]
                              [_ Variable:.:PI3Km -1]
                              [_ Variable:.:PIP3m 1]
                              [_ Variable:.:PI3Km 1];
      p 0.17;
    }



#identical dephosphorylation
  Process DiffusionInfluencedReactionProcess(PIP3toPIP2)
    {
      VariableReferenceList   [_ Variable:.:PIP3m -1]
                              [_ Variable:.:PTENm -1]
                              [_ Variable:.:PIP2m 1]
                              [_ Variable:.:PTENm 1];
      p 1;
    }
  Process DiffusionInfluencedReactionProcess(PIP3atoPIP2)
    {
      VariableReferenceList   [_ Variable:.:PIP3a -1]
                              [_ Variable:.:PTENm -1]
                              [_ Variable:.:PIP2m 1]
                              [_ Variable:.:PTENm 1];
      p 1;
    }
#identical dephosphorylation end


#Membrane dissociations:
  Process SpatiocyteNextReactionProcess(dissociatePTEN)
    {
      VariableReferenceList   [_ Variable:.:PTENm -1]
                              [_ Variable:.:PTEN 1];
      k 0.09;
    }
  Process SpatiocyteNextReactionProcess(dissociatePI3K)
    {
      VariableReferenceList   [_ Variable:.:PI3Km -1]
                              [_ Variable:.:PI3K 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP3)
    {
      VariableReferenceList   [_ Variable:.:PIP3m -1]
                              [_ Variable:.:PIP2 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP3a)
    {
      VariableReferenceList   [_ Variable:.:PIP3a -1]
                              [_ Variable:.:PIP2 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP2)
    {
      VariableReferenceList   [_ Variable:.:PIP2m -1]
                              [_ Variable:.:PIP2 1];
      k 0.0001;
    }
#Membrane dissociations end

  Process VisualizationLogProcess(logger)
    {
      VariableReferenceList [_ Variable:.:LIPID]
                            [_ Variable:.:PIP2m]
                            [_ Variable:.:PTENm]
                            [_ Variable:.:PIP3m]
                            [_ Variable:.:PIP3a]
                            [_ Variable:.:PI3Km];
                            
      LogInterval 20;
    }
  Process MoleculePopulateProcess(populate)
    {
      VariableReferenceList [_ Variable:.:PIP2m]
                            [_ Variable:.:PIP3m]
                            [_ Variable:.:PIP3a]
                            [_ Variable:.:PTENm]
                            [_ Variable:.:PI3Km];
    }
}




