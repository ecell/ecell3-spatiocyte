# A test model for Spatiocyte 
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 2e-8;
}

System System( / )
{
  StepperID       SS; 
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     1;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( LENGTHX )
    {
      Value     5e-6;        # in meters
    } 
  Variable Variable( LENGTHY )
    {
      Value     2e-6;        # in meters
    } 
  Variable Variable( LENGTHZ )
    {
      Value     5e-6;        # in meters
    } 
  Variable Variable( SIZE )
    {
      Value	 5e-14;  # volume in liters
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
      Value     0;               # { 0: Reflective (Impermeable) 
                                 #   1: Periodic (Permeable) }
    } 
  Variable Variable( SURFACEY )
    {
      Value     0;               # { 0: Reflective (Impermeable) 
                                 #   1: Periodic (Permeable) }
    } 
  Variable Variable( SURFACEZ )
    {
      Value     0;               # { 0: Reflective (Impermeable) 
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
  Variable Variable(A)
    {
      Value             1000; 
    } 
  Variable Variable(B)
    {
      Value             100; 
    } 
  Variable Variable(C)
    {
      Value             2000; 
    } 
  Process VisualizationLogProcess(loggerMean)
    {
      VariableReferenceList [_ Variable:/:A]
                            [_ Variable:/:B]
                            [_ Variable:/:C]
                            [_ Variable:/Membrane:D]
                            [_ Variable:/Membrane:E]
                            [_ Variable:/Membrane:F]
                            [_ Variable:/Membrane:LIPID];
      LogInterval 0.05;
      FileName "visualLog0.dat";
    }
}

System System( /Membrane )
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

  Variable Variable(D)
    {
      Value             1000; 
    } 
  Variable Variable(E)
    {
      Value             100; 
    } 
  Variable Variable(F)
    {
      Value             2000; 
    } 
}
