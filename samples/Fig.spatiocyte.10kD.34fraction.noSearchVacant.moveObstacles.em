# A test model for Spatiocyte 
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 3e-9;
  SearchVacant 0;
}

System System( / )
{
  StepperID       SS; 
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     2;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.5e-6;        # in meters
    } 
  Variable Variable( SIZE )
    {
      Value	 8.3333e-17;   # volume in liters
                              # RAM (memory) required volume/
                              # 0.85e-14 liter max for 8GB memory
                              # e.coli sphere radius: 0.5e-6 m 
                              # e.coli cylinder length: 3e-6 m
      # e.coli volume = (4*PI*(0.5e-6)^3)/3 + PI*(0.5e-6)^(2)*(3e-6)
      #               = 2.88e-18 m^3
      #               = 2.88e-15 l
    } 
  Variable Variable( SURFACEX )
    {
      Value     1;               # { 0: Reflective (Impermeable) 
                                 #   1: Periodic (Permeable) }
    } 
  Variable Variable( SURFACEY )
    {
      Value     1;               # { 0: Reflective (Impermeable) 
                                 #   1: Periodic (Permeable) }
    } 
  Variable Variable( SURFACEZ )
    {
      Value     1;               # { 0: Reflective (Impermeable) 
                                 #   1: Periodic (Permeable) }
    } 
  Variable Variable( ORIGINX )
    {
      Value     0;               
    } 
  Variable Variable( ORIGINY )
    {
      Value     0;          
    } 
  Variable Variable( ORIGINZ )
    {
      Value     0;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
  Variable Variable( E )
    {
      #200nM = Molar*(liter volume)*Avogadro
      #      = mol.l^{-1}.l.mol^{-1} 
      #      = (200e-9).(1e-15).(6.0221415e+23)
      #      = 120 molecules
      Value             1; 
    } 
  Variable Variable( S )
    {
      #1uM = Molar*(liter volume)*Avogadro
      #    = mol.l^{-1}.l.mol^{-1} 
      #    = (1e-6).(1e-15).(6.0221415e+23)
      #    = 602 molecules
      Value             5; 
    } 
  Variable Variable( ES )
    {
      Value             0; 
    } 
  Variable Variable( P )
    {
      Value             0; 
    } 
  Variable Variable( X )
    {
      Value             250512; 
    } 
  Process DiffusionProcess( diffuseX )
    {
      VariableReferenceList [ _ Variable:/:X ];
      D 0;
    }
  Process DiffusionProcess( diffuse )
    {
      VariableReferenceList [ _ Variable:/:ES ]
                            [ _ Variable:/:P ];
      D 2.5e-12;
    }
  Process DiffusionProcess( diffuse2 )
    {
      VariableReferenceList [ _ Variable:/:E ]
                            [ _ Variable:/:S ];
      D 2.5e-12;
    }
  Process DiffusionInfluencedReactionProcess( E_S_to_ES )
    {
      VariableReferenceList [_ Variable:/:E -1]
                            [_ Variable:/:S -1]
                            [_ Variable:/:ES 1];
      k 8e-18;
    }
  Process SpatiocyteNextReactionProcess( ES_to_E_P )
    {
      VariableReferenceList [_ Variable:/:ES -1]
                            [_ Variable:/:E 1]
                            [_ Variable:/:P 1];
      k 100;
    }
  Process SpatiocyteNextReactionProcess( ES_to_E_S )
    {
      VariableReferenceList [_ Variable:/:ES -1]
                            [_ Variable:/:E 1]
                            [_ Variable:/:S 1];
      k 100;
    }
  Process IteratingLogProcess( validate )
    {
      VariableReferenceList [_ Variable:/:E]
                            [_ Variable:/:S]
                            [_ Variable:/:ES]
                            [_ Variable:/:P];
      LogDuration 4;
      LogInterval 1e-3;
      Iterations 1000; 
      SaveInterval 1;
      FileName "Fig.spatiocyte.10kD.34fraction.noSearchVacant.csv";
    }
}



