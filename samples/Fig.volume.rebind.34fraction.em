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
      Value	 1e-15;   # volume in liters
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
      Value             1; 
    } 
  Variable Variable( X )
    {
      Value             3006260; 
    } 
  Process DiffusionProcess( diffuseX )
    {
      VariableReferenceList [ _ Variable:/:X ];
      D 0;
    }
  Process DiffusionProcess( diffuse )
    {
      VariableReferenceList [ _ Variable:/:E ]
                            [ _ Variable:/:S ];
      D 2.5e-12;
    }
  Process DiffusionInfluencedReactionProcess( E_S_to_ES )
    {
      VariableReferenceList [_ Variable:/:E -1]
                            [_ Variable:/:S -1]
                            [_ Variable:/:VACANT 1];
      k 8e-18;
    }
  Process IteratingLogProcess( validate )
    {
      VariableReferenceList [_ Variable:/:E];
      RebindTime 1;
      LogDuration 4;
      Iterations 100000000;
      SaveInterval 100;
      InContact 1;
      Centered 1;
      FileName "Fig.volume.rebind.34fraction.csv";
    }
}


