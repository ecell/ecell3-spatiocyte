# Wildtype MinDE model as published in
# Arjunan, S. N. V. and Tomita, M. (2010). A new multicompartmental
# reaction-diffusion modeling method links transient membrane attachment of
# E. coli MinE to E-ring formation. Syst. Synth. Biol. 4(1):35-53.
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 1e-8;    # m
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
      Value 1;         # { 0: Spherical (uses SIZE) 
                       #   1: Rod (uses SIZE, LENGTHY=2*radius)
                       #   2: Cubic (uses SIZE)
                       #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                       #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable(SIZE)
    {
      Value	3.27e-18;  # m^3
    } 
  Variable Variable(LENGTHY)
    {
      Value 1e-6;      # m
    } 
  Variable Variable(VACANT)
    {
      Value 0; 
    } 
  Variable Variable(MinDatp)
    {
      Value 0;         # molecule number 
    } 
  Variable Variable(MinDadp)
    {
      Value 1300;      # molecule number 
    } 
  Variable Variable(MinEE)
    {
      Value 0;         # molecule number 
    } 
  Process DiffusionProcess(diffuseMinD)
    {
      VariableReferenceList [_ Variable:/:MinDatp]
                            [_ Variable:/:MinDadp];
      D 16e-12;        # m^2/s
    }
  Process DiffusionProcess(diffuseMinE)
    {
      VariableReferenceList [_ Variable:/:MinEE];
      D 10e-12;        # m^2/s
    }
  Process VisualizationLogProcess(visualize)
    {
      VariableReferenceList [_ Variable:/Surface:MinEE]
                            [_ Variable:/Surface:MinDEE]
                            [_ Variable:/Surface:MinDEED]
                            [_ Variable:/Surface:MinD];
      LogInterval 0.5; # s
    }
  Process FluorescentImagingProcess(fluorescence)
    {
      VariableReferenceList [_ Variable:/Surface:MinEE 2]
                            [_ Variable:/Surface:MinDEE 3]
                            [_ Variable:/Surface:MinDEED 4]
                            [_ Variable:/Surface:MinD 1]
                            [_ Variable:/Surface:MinEE -2]
                            [_ Variable:/Surface:MinDEED -2]
                            [_ Variable:/Surface:MinEE -1]
                            [_ Variable:/Surface:MinDEED -4]
                            [_ Variable:/Surface:MinD -1];
      FileName "fluorescenceLog0.dat";
    }
  Process MoleculePopulateProcess(populate)
    {
      VariableReferenceList [_ Variable:/:MinDatp ]
                            [_ Variable:/:MinDadp ]
                            [_ Variable:/:MinEE ]
                            [_ Variable:/Surface:MinD ]
                            [_ Variable:/Surface:MinDEE ]
                            [_ Variable:/Surface:MinDEED ]
                            [_ Variable:/Surface:MinEE ];
    }
}

System System( /Surface )
{
  StepperID SS;

  Variable Variable( TYPE )
    {
      Value 1;         # { 0: Volume
                       #   1: Surface }
    } 
  Variable Variable( LIPID )
    {
      Value 0;
    } 
  Variable Variable(MinD)
    {
      Value 0;         # molecule number 
    } 
  Variable Variable(MinEE)
    {
      Value 0;         # molecule number 
    } 
  Variable Variable(MinDEE)
    {
      Value 700;       # molecule number 
    } 
  Variable Variable(MinDEED)
    {
      Value 0;         # molecule number 
    } 
  Process DiffusionProcess(diffuseMinD)
    {
      VariableReferenceList [_ Variable:/Surface:MinD];
      D 0.02e-12;      # m^2/s
    }
  Process DiffusionProcess(diffuseMinEE)
    {
      VariableReferenceList [_ Variable:/Surface:MinEE];
      D 0.02e-12;      # m^2/s
    }
  Process DiffusionProcess(diffuseMinDEE)
    {
      VariableReferenceList [_ Variable:/Surface:MinDEE];
      D 0.02e-12;      # m^2/s
    }
  Process DiffusionProcess(diffuseMinDEED)
    {
      VariableReferenceList [_ Variable:/Surface:MinDEED];
      D 0.02e-12;      # m^2/s
    }
  Process DiffusionInfluencedReactionProcess(reaction1) 
    {
      VariableReferenceList   [_ Variable:/Surface:LIPID -1]
                              [_ Variable:/:MinDatp -1]
                              [_ Variable:/Surface:MinD 1];
      k 2.2e-8;        # m/s
    } 
  Process DiffusionInfluencedReactionProcess(reaction2)
    {
      VariableReferenceList   [_ Variable:/Surface:MinD -1]
                              [_ Variable:/:MinDatp -1]
                              [_ Variable:/Surface:MinD 1]
                              [_ Variable:/Surface:MinD 1];
      k 3e-20;         # m^3/s
    }
  Process DiffusionInfluencedReactionProcess(reaction3)
    {
      VariableReferenceList   [_ Variable:/Surface:MinD -1]
                              [_ Variable:/:MinEE -1]
                              [_ Variable:/Surface:MinDEE 1];
      k 5e-19;         # m^3/s
    }
  Process SpatiocyteNextReactionProcess(reaction4)
    {
      VariableReferenceList   [_ Variable:/Surface:MinDEE -1]
                              [_ Variable:/Surface:MinEE 1]
                              [_ Variable:/:MinDadp 1];
      k 1;             # s^{-1}
    }
  Process SpatiocyteNextReactionProcess(reaction5)
    {
      VariableReferenceList [_ Variable:/:MinDadp -1]
                            [_ Variable:/:MinDatp 1];
      k 5;             # s^{-1}
    }
  Process DiffusionInfluencedReactionProcess(reaction6)
    {
      VariableReferenceList   [_ Variable:/Surface:MinDEE -1]
                              [_ Variable:/Surface:MinD -1]
                              [_ Variable:/Surface:MinDEED 1];
      k 5e-15;         # m^3/s
    }
  Process SpatiocyteNextReactionProcess(reaction7)
    {
      VariableReferenceList   [_ Variable:/Surface:MinDEED -1]
                              [_ Variable:/Surface:MinDEE 1]
                              [_ Variable:/:MinDadp 1];
      k 1;             # s^{-1}
    }
  Process SpatiocyteNextReactionProcess(reaction8)
    {
      VariableReferenceList   [_ Variable:/Surface:MinEE -1]
                              [_ Variable:/:MinEE 1];
      k 0.83;          # s^{-1}
    }
}

