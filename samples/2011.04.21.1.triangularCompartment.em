# Reaction between molecules from distinct surfaces
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

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
      Value 4e-5;      # m
    } 
  Variable Variable(LENGTHY)
    {
      Value 4e-5;      # m
    } 
  Variable Variable(LENGTHZ)
    {
      Value 0.3e-5;      # m
    } 
  Variable Variable(VACANT)
    {
      Value 0; 
    } 
  Process VisualizationLogProcess(loggerMean)
    {
      VariableReferenceList [_ Variable:/A/Cell/Surface:PIP2m]
                            [_ Variable:/A/Cell/Surface:PTENm]
                            [_ Variable:/A/Cell/Surface:PIP3m]
                            [_ Variable:/A/Cell/Surface:PIP3a]
                            [_ Variable:/A/Cell/Surface:PI3Km];
      LogInterval 1;
    }
  Process MoleculePopulateProcess(populate)
    {
      VariableReferenceList [_ Variable:/A/Cell/Surface:PIP2m]
                            [_ Variable:/A/Cell/Surface:PIP3m]
                            [_ Variable:/A/Cell/Surface:PIP3a]
                            [_ Variable:/A/Cell/Surface:PTENm]
                            [_ Variable:/A/Cell/Surface:PI3Km];
    }
}

System System(/A)
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
      Value 4e-5;      # m
    } 
  Variable Variable(LENGTHY)
    {
      Value 4e-5;      # m
    } 
  Variable Variable(LENGTHZ)
    {
      Value 0.22e-5;    # m
    } 
  Variable Variable(ORIGINZ)
    {
      Value 0;       
    } 
  Variable Variable(ORIGINY)
    {
      Value 0.7;       
    } 
  Variable Variable(ORIGINX)
    {
      Value 0.4;       
    } 
  Variable Variable(ROTATEZ)
    {
      Value 1.047;       
    } 
  Variable Variable(VACANT)
    {
      Value 0; 
    } 
}

System System(/A/Cell)
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
      Value 2.87e-5;      # m
    } 
  Variable Variable(LENGTHY)
    {
      Value 2.49e-5;      # m
    } 
  Variable Variable(LENGTHZ)
    {
      Value 0.1e-5;    # m
    } 
  Variable Variable(ROTATEZ)
    {
      Value 2.09439;       
    } 
  Variable Variable(ORIGINX)
    {
      Value -0.9;       
    } 
  Variable Variable(ORIGINY)
    {
      Value 0;       
    } 
  Variable Variable(VACANT)
    {
      Value 0; 
    } 
}

System System(/A/Surface)
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
}

System System(/A/Cell/Surface)
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
  Variable Variable(PIP2m)
    {
      Value             1000; 
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
      Value             250; 
    } 
  Variable Variable(PI3Km)
    {
      Value             2500; 
    } 
  Variable Variable(PIP2)
    {
      Value             0.65e+4; 
      Name "HD";
    } 
  Variable Variable(PI3K)
    {
      Value             0.75e+4; 
      Name "HD";
    } 
  Variable Variable(PTEN)
    {
      Value             5000; 
      Name "HD";
    }
  Process DiffusionProcess(diffusePIP2)
    {
      VariableReferenceList [_ Variable:/A/Cell/Surface:PIP2m];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePIP3)
    {
      VariableReferenceList [_ Variable:/A/Cell/Surface:PIP3m];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePIP3a)
    {
      VariableReferenceList [_ Variable:/A/Cell/Surface:PIP3a];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePTEN)
    {
      VariableReferenceList [_ Variable:/A/Cell/Surface:PTENm];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePI3K)
    {
      VariableReferenceList [_ Variable:/A/Cell/Surface:PI3Km];
      D 1e-14;
    }

#Membrane recruitments:
  Process SpatiocyteNextReactionProcess(recruitPIP2)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PIP2 -1]
                              [_ Variable:/A/Cell/Surface:PIP2m 1];
      k 4e-2;
    }
  Process SpatiocyteNextReactionProcess(recruitPTEN)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PTEN -1]
                              [_ Variable:/A/Cell/Surface:PIP2m -1]
                              [_ Variable:/A/Cell/Surface:PTENm 1]
                              [_ Variable:/A/Cell/Surface:PIP2m 1];
      k 2e-14;
    }
  Process SpatiocyteNextReactionProcess(recruitPI3Ka)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PIP3a -1]
                              [_ Variable:/A/Cell/Surface:PI3K -1]
                              [_ Variable:/A/Cell/Surface:PIP3m 1]
                              [_ Variable:/A/Cell/Surface:PI3Km 1];
      k 1e-13;
    }
#Membrane recruitments end


#Activations:
  Process DiffusionInfluencedReactionProcess(dimerPIP3)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PIP3m -1]
                              [_ Variable:/A/Cell/Surface:PIP3m -1]
                              [_ Variable:/A/Cell/Surface:PIP3a 1]
                              [_ Variable:/A/Cell/Surface:PIP3a 1];
      p 0.65;
    }
  Process DiffusionInfluencedReactionProcess(PIP2toPIP3)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PIP2m -1]
                              [_ Variable:/A/Cell/Surface:PI3Km -1]
                              [_ Variable:/A/Cell/Surface:PIP3m 1]
                              [_ Variable:/A/Cell/Surface:PI3Km 1];
      p 0.17;
    }



#identical dephosphorylation
  Process DiffusionInfluencedReactionProcess(PIP3toPIP2)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PIP3m -1]
                              [_ Variable:/A/Cell/Surface:PTENm -1]
                              [_ Variable:/A/Cell/Surface:PIP2m 1]
                              [_ Variable:/A/Cell/Surface:PTENm 1];
      p 1;
    }
  Process DiffusionInfluencedReactionProcess(PIP3atoPIP2)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PIP3a -1]
                              [_ Variable:/A/Cell/Surface:PTENm -1]
                              [_ Variable:/A/Cell/Surface:PIP2m 1]
                              [_ Variable:/A/Cell/Surface:PTENm 1];
      p 1;
    }
#identical dephosphorylation end


#Membrane dissociations:
  Process SpatiocyteNextReactionProcess(dissociatePTEN)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PTENm -1]
                              [_ Variable:/A/Cell/Surface:PTEN 1];
      k 0.09;
    }
  Process SpatiocyteNextReactionProcess(dissociatePI3K)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PI3Km -1]
                              [_ Variable:/A/Cell/Surface:PI3K 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP3)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PIP3m -1]
                              [_ Variable:/A/Cell/Surface:PIP2 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP3a)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PIP3a -1]
                              [_ Variable:/A/Cell/Surface:PIP2 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP2)
    {
      VariableReferenceList   [_ Variable:/A/Cell/Surface:PIP2m -1]
                              [_ Variable:/A/Cell/Surface:PIP2 1];
      k 0.0001;
    }
#Membrane dissociations end

}


