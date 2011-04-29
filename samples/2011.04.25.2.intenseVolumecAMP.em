Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 6e-8;
}

System System( / )
{
  StepperID       SS; 
  Variable Variable( SHAPE )
    {
      Value     3;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( LENGTHX )
    {
      Value     25e-6;        # in meters
    } 
  Variable Variable( LENGTHY )
    {
      Value     25e-6;        # in meters
    } 
  Variable Variable( LENGTHZ )
    {
      Value     1e-6;        # in meters
    } 
  Variable Variable( YZPLANE )
    {
      Value     5;               # { 0: Reflective (Impermeable) 
                                 #   1: Periodic (Permeable) }
    } 
  Variable Variable( XZPLANE )
    {
      Value     5;               # { 0: Reflective (Impermeable) 
                                 #   1: Periodic (Permeable) }
    } 
  Variable Variable( XYPLANE )
    {
      Value     3;               # { 0: Reflective (Impermeable) 
                                 #   1: Periodic (Permeable) }
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
  Variable Variable( cAMP )
    {
      Value             0; 
    } 
  Process DiffusionProcess(diffusecAMP)
    {
      VariableReferenceList [_ Variable:/:cAMP];
      D 1e-14;
    }
  Process DiffusionProcess(diffusecAMPblock)
    {
      VariableReferenceList [_ Variable:/Block:cAMP];
      D 1e-14;
    }
  Process SpatiocyteNextReactionProcess(degradecAMP)
    {
      VariableReferenceList   [_ Variable:/:cAMP -1]
                              [_ Variable:/:VACANT 1];
      k 0.04;
    }
  Process SpatiocyteNextReactionProcess(degradecAMPblock)
    {
      VariableReferenceList   [_ Variable:/Block:cAMP -1]
                              [_ Variable:/Block:VACANT 1];
      k 0.04;
    }
  Process DiffusionInfluencedReactionProcess(recruitcAMP)
    {
      VariableReferenceList [_ Variable:/:cAMP -1]
                            [_ Variable:/Membrane:LIPID -1]
                            [_ Variable:/Membrane:cAMP 1];
      p 0.6;
    }
  Process VisualizationLogProcess(loggerMean)
    {
      VariableReferenceList [_ Variable:/Membrane:PIP2m]
                            [_ Variable:/Membrane:PTENm]
                            [_ Variable:/Membrane:PIP3m]
                            [_ Variable:/Membrane:PIP3a]
                            [_ Variable:/Membrane:PI3Km]
                            [_ Variable:/Block:cAMP]
                            [_ Variable:/Membrane:cAMP]
                            [_ Variable:/:cAMP];
      LogInterval 10;
    }
  Process MoleculePopulateProcess( populate )
    {
      VariableReferenceList [_ Variable:/Membrane:PIP2m]
                            [_ Variable:/Membrane:PIP3m]
                            [_ Variable:/Membrane:PIP3a]
                            [_ Variable:/Membrane:PTENm]
                            [_ Variable:/Membrane:PI3Km]
                            [_ Variable:/Block:cAMP]
                            [_ Variable:/Membrane:cAMP]
                            [_ Variable:/:cAMP];
    }
}

System System( /Block )
{
  StepperID SS;
  Variable Variable( SHAPE )
    {
      Value     3;     
    } 
  Variable Variable( LENGTHX )
    {
      Value     1e-6;        # in meters
    } 
  Variable Variable( LENGTHY )
    {
      Value     26e-6;        # in meters
    } 
  Variable Variable( LENGTHZ )
    {
      Value     26e-6;        # in meters
    } 
  Variable Variable( ORIGINY )
    {
      Value     5; 
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
  Variable Variable( cAMP )
    {
      Value             0; 
    } 
  Process DiffusionInfluencedReactionProcess(interCompartmentDiffuse)
    {
      VariableReferenceList [_ Variable:/:cAMP -1]
                            [_ Variable:/Block:VACANT -1]
                            [_ Variable:/Block:cAMP 1];
      p 1;
    }
  Process DiffusionInfluencedReactionProcess(phosphataseBindcAMP)
    {
      VariableReferenceList [_ Variable:/Block:cAMP -1]
                            [_ Variable:/:VACANT -1]
                            [_ Variable:/:cAMP 1];
      p 1;
    }
}


System System( /Membrane )
{
  StepperID SS;

  Variable Variable(TYPE)
    {
      Value 1;               # { 0: Cytoplasm (requires SHAPE)
                             #   1: Surface }
    } 
  Variable Variable(LIPID)
    {
      Value 0;
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
  Variable Variable(cAMP)
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
      VariableReferenceList [_ Variable:/Membrane:PIP2m];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePIP3)
    {
      VariableReferenceList [_ Variable:/Membrane:PIP3m];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePIP3a)
    {
      VariableReferenceList [_ Variable:/Membrane:PIP3a];
      D 1e-14;
    }
  Process DiffusionProcess(diffusecAMP)
    {
      VariableReferenceList [_ Variable:/Membrane:cAMP];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePTEN)
    {
      VariableReferenceList [_ Variable:/Membrane:PTENm];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePI3K)
    {
      VariableReferenceList [_ Variable:/Membrane:PI3Km];
      D 1e-14;
    }


#Membrane recruitments:
  Process SpatiocyteNextReactionProcess(recruitPIP2)
    {
      VariableReferenceList   [_ Variable:/Membrane:PIP2 -1]
                              [_ Variable:/Membrane:PIP2m 1];
      k 4e-2;
    }
  Process SpatiocyteNextReactionProcess(recruitPTEN)
    {
      VariableReferenceList   [_ Variable:/Membrane:PTEN -1]
                              [_ Variable:/Membrane:PIP2m -1]
                              [_ Variable:/Membrane:PTENm 1]
                              [_ Variable:/Membrane:PIP2m 1];
      k 2e-14;
    }
  Process SpatiocyteNextReactionProcess(recruitPI3Ka)
    {
      VariableReferenceList   [_ Variable:/Membrane:PIP3a -1]
                              [_ Variable:/Membrane:PI3K -1]
                              [_ Variable:/Membrane:PIP3m 1]
                              [_ Variable:/Membrane:PI3Km 1];
      k 1e-13;
    }
  Process SpatiocyteNextReactionProcess(recruitPI3Kc)
    {
      VariableReferenceList   [_ Variable:/Membrane:cAMP -1]
                              [_ Variable:/Membrane:PI3K -1]
                              [_ Variable:/Membrane:PI3Km 1];
      k 1e-17;
    }
#Membrane recruitments end


#Activations:
  Process DiffusionInfluencedReactionProcess(dimerPIP3)
    {
      VariableReferenceList   [_ Variable:/Membrane:PIP3m -1]
                              [_ Variable:/Membrane:PIP3m -1]
                              [_ Variable:/Membrane:PIP3a 1]
                              [_ Variable:/Membrane:PIP3a 1];
      p 0.65;
    }
  Process DiffusionInfluencedReactionProcess(PIP2toPIP3)
    {
      VariableReferenceList   [_ Variable:/Membrane:PIP2m -1]
                              [_ Variable:/Membrane:PI3Km -1]
                              [_ Variable:/Membrane:PIP3m 1]
                              [_ Variable:/Membrane:PI3Km 1];
      p 0.17;
    }



#identical dephosphorylation
  Process DiffusionInfluencedReactionProcess(PIP3toPIP2)
    {
      VariableReferenceList   [_ Variable:/Membrane:PIP3m -1]
                              [_ Variable:/Membrane:PTENm -1]
                              [_ Variable:/Membrane:PIP2m 1]
                              [_ Variable:/Membrane:PTENm 1];
      p 1;
    }
  Process DiffusionInfluencedReactionProcess(PIP3atoPIP2)
    {
      VariableReferenceList   [_ Variable:/Membrane:PIP3a -1]
                              [_ Variable:/Membrane:PTENm -1]
                              [_ Variable:/Membrane:PIP2m 1]
                              [_ Variable:/Membrane:PTENm 1];
      p 1;
    }
#identical dephosphorylation end


#Membrane dissociations:
  Process SpatiocyteNextReactionProcess(dissociatePTEN)
    {
      VariableReferenceList   [_ Variable:/Membrane:PTENm -1]
                              [_ Variable:/Membrane:PTEN 1];
      k 0.09;
    }
  Process SpatiocyteNextReactionProcess(dissociatePI3K)
    {
      VariableReferenceList   [_ Variable:/Membrane:PI3Km -1]
                              [_ Variable:/Membrane:PI3K 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP3)
    {
      VariableReferenceList   [_ Variable:/Membrane:PIP3m -1]
                              [_ Variable:/Membrane:PIP2 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP3a)
    {
      VariableReferenceList   [_ Variable:/Membrane:PIP3a -1]
                              [_ Variable:/Membrane:PIP2 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP2)
    {
      VariableReferenceList   [_ Variable:/Membrane:PIP2m -1]
                              [_ Variable:/Membrane:PIP2 1];
      k 0.0001;
    }
#Membrane dissociations end
    Process SpatiocyteNextReactionProcess(AreleasecAMP)
    {
      VariableReferenceList   [_ Variable:/Membrane:PIP3m -1]
                              [_ Variable:/Membrane:PIP3m 1]
                              [_ Variable:/:cAMP 1];
      k 1;
    }
  Process SpatiocyteNextReactionProcess(dissociatecAMP)
    {
      VariableReferenceList   [_ Variable:/Membrane:cAMP -1]
                              [_ Variable:/Membrane:LIPID 1];
      k 0.5;
    }
}


