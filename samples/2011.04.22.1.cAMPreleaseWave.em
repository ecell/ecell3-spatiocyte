Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 6e-8;
}

System System( / )
{
  StepperID       SS; 
  Variable Variable( SHAPE )
    {
      Value     3;   
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
      Value     6e-6;        # in meters
    } 
  Variable Variable( YZPLANE )
    {
      Value     5;  
    } 
  Variable Variable( XZPLANE )
    {
      Value     5; 
    } 
  Variable Variable( XYPLANE )
    {
      Value     3;
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
  Process VisualizationLogProcess(loggerMean)
    {
      VariableReferenceList [_ Variable:/CellA/Surface:PIP2m]
                            [_ Variable:/CellA/Surface:PTENm]
                            [_ Variable:/CellA/Surface:PIP3m]
                            [_ Variable:/CellA/Surface:PIP3a]
                            [_ Variable:/CellA/Surface:PI3Km]
                            [_ Variable:/Membrane:cAMP];
      LogInterval 20;
    }
  Process MoleculePopulateProcess( populate )
    {
      VariableReferenceList [_ Variable:/CellA/Surface:PIP2m]
                            [_ Variable:/CellA/Surface:PTENm]
                            [_ Variable:/CellA/Surface:PIP3m]
                            [_ Variable:/CellA/Surface:PIP3a]
                            [_ Variable:/CellA/Surface:PI3Km]
                            [_ Variable:/Membrane:cAMP];
    }
}

System System( /Membrane )
{
  StepperID SS;

  Variable Variable(TYPE)
    {
      Value 1; 
    } 
  Variable Variable(LIPID)
    {
      Value 0;
    } 
  Variable Variable(cAMP)
    {
      Value 0;
    } 
  Variable Variable(REACTIVE)
    {
      Name "/CellA:Surface";
    }
  Process DiffusionProcess(diffusecAMP)
    {
      VariableReferenceList [_ Variable:/Membrane:cAMP];
      D 0.2e-12;
    }
  Process SpatiocyteNextReactionProcess(degradecAMP)
    {
      VariableReferenceList   [_ Variable:/Membrane:cAMP -1]
                              [_ Variable:/Membrane:LIPID 1];
      k 0.01;
    }
}

System System( /CellA )
{
  StepperID SS;
  Variable Variable( SHAPE )
    {
      Value     4;     
    } 
  Variable Variable( LENGTHX )
    {
      Value     15e-6;        # in meters
    } 
  Variable Variable( LENGTHY )
    {
      Value     15e-6;        # in meters
    } 
  Variable Variable( LENGTHZ )
    {
      Value     6e-6;        # in meters
    } 
  Variable Variable( ORIGINZ )
    {
      Value     -0.7;       
    } 
  Variable Variable( ORIGINX )
    {
      Value     0;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
}

System System( /CellA/Surface )
{
  StepperID SS;

  Variable Variable( TYPE )
    {
      Value 1;     
    } 
  Variable Variable( LIPID )
    {
      Value 1; # 1: enclosed surface
    } 
  Variable Variable(PIP2m)
    {
      Value             686; 
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
      Value             172; 
    } 
  Variable Variable(PI3Km)
    {
      Value             1722; 
    } 
  Variable Variable(PIP2)
    {
      Value             0.4478e+4; 
      Name "HD";
    } 
  Variable Variable(PI3K)
    {
      Value             0.515e+4; 
      Name "HD";
    } 
  Variable Variable(PTEN)
    {
      Value             3444; 
      Name "HD";
    } 
  Process DiffusionProcess(AdiffusePIP2)
    {
      VariableReferenceList [_ Variable:/CellA/Surface:PIP2m];
      D 1e-14;
    }
  Process DiffusionProcess(AdiffusePIP3)
    {
      VariableReferenceList [_ Variable:/CellA/Surface:PIP3m];
      D 1e-14;
    }
  Process DiffusionProcess(AdiffusePIP3a)
    {
      VariableReferenceList [_ Variable:/CellA/Surface:PIP3a];
      D 1e-14;
    }
  Process DiffusionProcess(AdiffusePTEN)
    {
      VariableReferenceList [_ Variable:/CellA/Surface:PTENm];
      D 1e-14;
    }
  Process DiffusionProcess(AdiffusePI3K)
    {
      VariableReferenceList [_ Variable:/CellA/Surface:PI3Km];
      D 1e-14;
    }

#Membrane recruitments:
  Process SpatiocyteNextReactionProcess(ArecruitPIP2)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PIP2 -1]
                              [_ Variable:/CellA/Surface:PIP2m 1];
      k 4e-2;
    }
  Process SpatiocyteNextReactionProcess(ArecruitPTEN)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PTEN -1]
                              [_ Variable:/CellA/Surface:PIP2m -1]
                              [_ Variable:/CellA/Surface:PTENm 1]
                              [_ Variable:/CellA/Surface:PIP2m 1];
      k 2e-14;
    }
  Process SpatiocyteNextReactionProcess(ArecruitPI3Ka)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PIP3a -1]
                              [_ Variable:/CellA/Surface:PI3K -1]
                              [_ Variable:/CellA/Surface:PIP3m 1]
                              [_ Variable:/CellA/Surface:PI3Km 1];
      k 1e-13;
    }
#Membrane recruitments end


#Activations:
  Process DiffusionInfluencedReactionProcess(AdimerPIP3)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PIP3m -1]
                              [_ Variable:/CellA/Surface:PIP3m -1]
                              [_ Variable:/CellA/Surface:PIP3a 1]
                              [_ Variable:/CellA/Surface:PIP3a 1];
      p 0.65;
    }
  Process DiffusionInfluencedReactionProcess(APIP2toPIP3)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PIP2m -1]
                              [_ Variable:/CellA/Surface:PI3Km -1]
                              [_ Variable:/CellA/Surface:PIP3m 1]
                              [_ Variable:/CellA/Surface:PI3Km 1];
      p 0.17;
    }



#identical dephosphorylation
  Process DiffusionInfluencedReactionProcess(APIP3toPIP2)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PIP3m -1]
                              [_ Variable:/CellA/Surface:PTENm -1]
                              [_ Variable:/CellA/Surface:PIP2m 1]
                              [_ Variable:/CellA/Surface:PTENm 1];
      p 1;
    }
  Process DiffusionInfluencedReactionProcess(APIP3atoPIP2)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PIP3a -1]
                              [_ Variable:/CellA/Surface:PTENm -1]
                              [_ Variable:/CellA/Surface:PIP2m 1]
                              [_ Variable:/CellA/Surface:PTENm 1];
      p 1;
    }
#identical dephosphorylation end


#Membrane dissociations:
  Process SpatiocyteNextReactionProcess(AdissociatePTEN)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PTENm -1]
                              [_ Variable:/CellA/Surface:PTEN 1];
      k 0.09;
    }
  Process SpatiocyteNextReactionProcess(AdissociatePI3K)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PI3Km -1]
                              [_ Variable:/CellA/Surface:PI3K 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(AdissociatePIP3)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PIP3m -1]
                              [_ Variable:/CellA/Surface:PIP2 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(AdissociatePIP3a)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PIP3a -1]
                              [_ Variable:/CellA/Surface:PIP2 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(AdissociatePIP2)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PIP2m -1]
                              [_ Variable:/CellA/Surface:PIP2 1];
      k 0.0001;
    }
#Membrane dissociations end

  Process SpatiocyteNextReactionProcess(AreleasecAMP)
    {
      VariableReferenceList   [_ Variable:/CellA/Surface:PIP3m -1]
                              [_ Variable:/CellA/Surface:PIP3m 1]
                              [_ Variable:/Membrane:cAMP 1];
      k 0.09;
    }
}

