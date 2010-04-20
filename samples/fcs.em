# A test model for Spatiocyte 
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 2e-9;
  SearchVacant 1;
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
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ)
                                #   3: Cuboid (uses LENGTHY, LENGTHZ) }
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
      Value             500; 
    } 
  Variable Variable( X )
    {
      Value             10146127; 
    } 
  Process DiffusionProcess( diffuseX )
    {
      VariableReferenceList [ _ Variable:/:X ];
      D 0;
    } 
  Process DiffusionProcess( diffuse )
    {
      VariableReferenceList [ _ Variable:/:E ];
      D 24e-12;
    }
  Process CoordinateLogProcess( coordinate )
    {
      VariableReferenceList [_ Variable:/:E];
      FileName "coordinates.csv";
      LogDuration 100;
    }
}


