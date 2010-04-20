# A test model for Spatiocyte 
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 10e-9;
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
      Value     1;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( LENGTHX )
    {
      Value     20e-6;        # in meters
    } 
  Variable Variable( LENGTHY )
    {
      Value     1e-6;        # in meters
    } 
  Variable Variable( LENGTHZ )
    {
      Value     20e-6;        # in meters
    } 
  Variable Variable( SIZE )
    {
      Value	 3.27380952e-15;   # volume in liters
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
  Process VisualizationLogProcess(logger)
    {
      VariableReferenceList [_ Variable:/Surface:A ]
                            [_ Variable:/Surface:B ]
                            [_ Variable:/Surface:C ]
                            [_ Variable:/Surface:D ];
      LogInterval 0.05;
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
  Variable Variable( A )
    {
      Value             0; 
    } 
  Variable Variable( B )
    {
      Value             20; 
    } 
  Variable Variable( C )
    {
      Value             5000; 
    } 
  Variable Variable( D )
    {
      Value             0; 
    } 
  Process DiffusionProcess(diffuse)
    {
      VariableReferenceList [_ Variable:/Surface:A]
                            [_ Variable:/Surface:C];
      D 2.5e-12;
    }
  Process PolymerizationParameterProcess(poly)
    {
      VariableReferenceList [_ Variable:/Surface:B]
                            [_ Variable:/Surface:D];
      BendAngles [0.77] [0];
      PolymerDirectionality 0;
    }
  Process PolymerizationProcess(poly3)
    {
      VariableReferenceList [_ Variable:/Surface:A -1]
                            [_ Variable:/Surface:A -1]
                            [_ Variable:/Surface:B 1]
                            [_ Variable:/Surface:B 1];
      BendAngle 0.77;
      p 1;
    }
  Process PolymerizationProcess(poly4)
    {
      VariableReferenceList [_ Variable:/Surface:B -1]
                            [_ Variable:/Surface:C -1]
                            [_ Variable:/Surface:D 1]
                            [_ Variable:/Surface:B 1];
      BendAngle 0.77;
      p 1;
    }
}





