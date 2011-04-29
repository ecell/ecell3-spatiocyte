# A test model for Spatiocyte 
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 6e-8;
}

System System( / )
{
  StepperID       SS; 
  Variable Variable( SHAPE )
    {
      Value     4;              #{ 0: Spherical
                                #  1: Rod (uses RADIUS)
                                #  2: Cubic (uses SurfaceX,SurfaceY,SurfaceZ)
                                #  3: Cuboid (uses SurfaceX,SurfaceY,SurfaceZ) }
    } 
  Variable Variable( LENGTHX )
    {
      Value     20e-6;        # in meters
    } 
  Variable Variable( LENGTHY )
    {
      Value     20e-6;        # in meters
    } 
  Variable Variable( LENGTHZ )
    {
      Value     20e-6;        # in meters
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
  Process VisualizationLogProcess(loggerMean)
    {
      VariableReferenceList [_ Variable:/Surface:PIP2m]
                            [_ Variable:/Surface:PTENm]
                            [_ Variable:/Surface:PIP3m]
                            [_ Variable:/Surface:PIP3a]
                            [_ Variable:/Surface:PI3Km];
      LogInterval 20;
      FileName "visualLogSphereHD0.dat";
    }
  Process MoleculePopulateProcess( populate )
    {
      VariableReferenceList [_ Variable:/Surface:PIP2m]
                            [_ Variable:/Surface:PIP3m]
                            [_ Variable:/Surface:PIP3a]
                            [_ Variable:/Surface:PTENm]
                            [_ Variable:/Surface:PI3Km];
    }
}

System System( /Surface )
{
  StepperID SS;

  Variable Variable( TYPE )
    {
      Value 1;               # { 0: Cytoplasm (requires SHAPE)
                             #   1: Surface }
    } 
  Variable Variable( LIPID )
    {
      Value 0; # value is always 0 for lipids;
    } 
  Variable Variable(PIP2m)
    {
      Value             1840; 
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
      Value             461; 
    } 
  Variable Variable(PI3Km)
    {
      Value             4614; 
    } 
  Variable Variable(PIP2)
    {
      Value             1.2e+4; 
      Name "HD";
    } 
  Variable Variable(PI3K)
    {
      Value             1.38e+4; 
      Name "HD";
    } 
  Variable Variable(PTEN)
    {
      Value             9229; 
      Name "HD";
    } 
  Process DiffusionProcess(diffusePIP2)
    {
      VariableReferenceList [_ Variable:/Surface:PIP2m];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePIP3)
    {
      VariableReferenceList [_ Variable://Surface:PIP3m];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePIP3a)
    {
      VariableReferenceList [_ Variable:/Surface:PIP3a];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePTEN)
    {
      VariableReferenceList [_ Variable:/Surface:PTENm];
      D 1e-14;
    }
  Process DiffusionProcess(diffusePI3K)
    {
      VariableReferenceList [_ Variable:/Surface:PI3Km];
      D 1e-14;
    }

#Membrane recruitments:
  Process SpatiocyteNextReactionProcess(recruitPIP2)
    {
      VariableReferenceList   [_ Variable:/Surface:PIP2 -1]
                              [_ Variable:/Surface:PIP2m 1];
      k 4e-2;
    }
  Process SpatiocyteNextReactionProcess(recruitPTEN)
    {
      VariableReferenceList   [_ Variable:/Surface:PTEN -1]
                              [_ Variable:/Surface:PIP2m -1]
                              [_ Variable:/Surface:PTENm 1]
                              [_ Variable:/Surface:PIP2m 1];
      k 2e-14;
    }
  Process SpatiocyteNextReactionProcess(recruitPI3Ka)
    {
      VariableReferenceList   [_ Variable:/Surface:PIP3a -1]
                              [_ Variable:/Surface:PI3K -1]
                              [_ Variable:/Surface:PIP3m 1]
                              [_ Variable:/Surface:PI3Km 1];
      k 1e-13;
    }
#Membrane recruitments end


#Activations:
  Process DiffusionInfluencedReactionProcess(dimerPIP3)
    {
      VariableReferenceList   [_ Variable:/Surface:PIP3m -1]
                              [_ Variable:/Surface:PIP3m -1]
                              [_ Variable:/Surface:PIP3a 1]
                              [_ Variable:/Surface:PIP3a 1];
      p 0.65;
    }
  Process DiffusionInfluencedReactionProcess(PIP2toPIP3)
    {
      VariableReferenceList   [_ Variable:/Surface:PIP2m -1]
                              [_ Variable:/Surface:PI3Km -1]
                              [_ Variable:/Surface:PIP3m 1]
                              [_ Variable:/Surface:PI3Km 1];
      p 0.17;
    }



#identical dephosphorylation
  Process DiffusionInfluencedReactionProcess(PIP3toPIP2)
    {
      VariableReferenceList   [_ Variable:/Surface:PIP3m -1]
                              [_ Variable:/Surface:PTENm -1]
                              [_ Variable:/Surface:PIP2m 1]
                              [_ Variable:/Surface:PTENm 1];
      p 1;
    }
  Process DiffusionInfluencedReactionProcess(PIP3atoPIP2)
    {
      VariableReferenceList   [_ Variable:/Surface:PIP3a -1]
                              [_ Variable:/Surface:PTENm -1]
                              [_ Variable:/Surface:PIP2m 1]
                              [_ Variable:/Surface:PTENm 1];
      p 1;
    }
#identical dephosphorylation end


#Membrane dissociations:
  Process SpatiocyteNextReactionProcess(dissociatePTEN)
    {
      VariableReferenceList   [_ Variable:/Surface:PTENm -1]
                              [_ Variable:/Surface:PTEN 1];
      k 0.09;
    }
  Process SpatiocyteNextReactionProcess(dissociatePI3K)
    {
      VariableReferenceList   [_ Variable:/Surface:PI3Km -1]
                              [_ Variable:/Surface:PI3K 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP3)
    {
      VariableReferenceList   [_ Variable:/Surface:PIP3m -1]
                              [_ Variable:/Surface:PIP2 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP3a)
    {
      VariableReferenceList   [_ Variable:/Surface:PIP3a -1]
                              [_ Variable:/Surface:PIP2 1];
      k 0.02;
    }
  Process SpatiocyteNextReactionProcess(dissociatePIP2)
    {
      VariableReferenceList   [_ Variable:/Surface:PIP2m -1]
                              [_ Variable:/Surface:PIP2 1];
      k 0.0001;
    }
#Membrane dissociations end

}




