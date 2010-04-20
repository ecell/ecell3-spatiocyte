# A test model for Spatiocyte 
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 2.5e-9;
}

System System( / )
{
  StepperID       SS; 
  Variable Variable( TYPE )
    {
      Value     0;               # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                 #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     2;               # { 0: Spherical
                                 #   1: Rod (uses RADIUS)
                                 #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.1e-6;          # in meters
    } 
  Variable Variable( SIZE )
    {
      Value	 1e-15;    # volume in liters
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
      Value     0; 
    } 
  Variable Variable( A )
    {
      Value             1; 
    } 
  Variable Variable( B )
    {
      Value             1; 
    } 
  Process DiffusionProcess( diffuseD0_1 )
    {
      VariableReferenceList [ _ Variable:/:A ]
                            [ _ Variable:/:B ];
      D 1e-12;
    }
  Process DiffusionInfluencedReactionProcess( EF )
    {
      VariableReferenceList [ _ Variable:/:A -1]
                            [ _ Variable:/:B -1]
                            [ _ Variable:/:VACANT 1];

      k 0.092e-18;
    }
  Process IteratingLogProcess( validate )
    {
      VariableReferenceList [_ Variable:/:A];
      RebindTime 1;
      LogDuration 100;
      Iterations 100000000;
      SaveInterval 100;
      InContact 1;
      Centered 1;
      FileName "Fig.volume.irreversible.csv";
    }
}


