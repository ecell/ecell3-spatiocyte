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
      Value     4;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( LENGTHX )
    {
      Value     8e-6;        # in meters
    } 
  Variable Variable( LENGTHY )
    {
      Value     3.5e-6;        # in meters
    } 
  Variable Variable( LENGTHZ )
    {
      Value     8e-6;        # in meters
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
                            [_ Variable:/Membrane:LIPID]
                            [_ Variable:/Nucleus:G]
                            [_ Variable:/Nucleus:H]
                            [_ Variable:/Nucleus:I]
                            [_ Variable:/Nucleus/Membrane:J]
                            [_ Variable:/Nucleus/Membrane:K]
                            [_ Variable:/Nucleus/Membrane:L]
                            [_ Variable:/Nucleus/Membrane:LIPID]
                            [_ Variable:/MitochondrionA/Membrane:LIPID]
                            [_ Variable:/MitochondrionB/Membrane:LIPID]
                            [_ Variable:/MitochondrionC/Membrane:LIPID]
                            [_ Variable:/Nucleus/Nucleolus:DNA]
                            [_ Variable:/VacuoleA/Membrane:LIPID]
                            [_ Variable:/VacuoleB/Membrane:LIPID]
                            [_ Variable:/VacuoleC/Membrane:LIPID]
                            [_ Variable:/LysosomeA/Membrane:LIPID]
                            [_ Variable:/LysosomeB/Membrane:LIPID]
                            [_ Variable:/LysosomeC/Membrane:LIPID];
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

System System( /Nucleus )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( SIZE )
    {
      Value	 5e-15;  # volume in liters
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
  Variable Variable(G)
    {
      Value             1000; 
    } 
  Variable Variable(H)
    {
      Value             100; 
    } 
  Variable Variable(I)
    {
      Value             2000; 
    } 
}

System System( /MitochondrionA )
{
  StepperID SS;
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
  Variable Variable( LENGTHY )
    {
      Value     0.4e-6;        # in meters
    } 
  Variable Variable( SIZE )
    {
      Value	 0.15e-15;  # volume in liters
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
      Value     0.60;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
}

System System( /MitochondrionA/Membrane )
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
}

System System( /MitochondrionB )
{
  StepperID SS;
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
  Variable Variable( LENGTHY )
    {
      Value     0.4e-6;        # in meters
    } 
  Variable Variable( SIZE )
    {
      Value	 0.15e-15;  # volume in liters
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
      Value     -0.3;          
    } 
  Variable Variable( ORIGINZ )
    {
      Value     -0.60;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
}

System System( /MitochondrionB/Membrane )
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
}

System System( /MitochondrionC )
{
  StepperID SS;
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
  Variable Variable( LENGTHY )
    {
      Value     0.4e-6;        # in meters
    } 
  Variable Variable( SIZE )
    {
      Value	 0.15e-15;  # volume in liters
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
      Value     0.2;               
    } 
  Variable Variable( ORIGINY )
    {
      Value     0.6;          
    } 
  Variable Variable( ORIGINZ )
    {
      Value     -0.2;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
}

System System( /MitochondrionC/Membrane )
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
}

System System( /Nucleus/Membrane )
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
  Variable Variable(J)
    {
      Value             1000; 
    } 
  Variable Variable(K)
    {
      Value             100; 
    } 
  Variable Variable(L)
    {
      Value             2000; 
    } 
}


System System( /Nucleus/Nucleolus )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( SIZE )
    {
      Value	 0.15e-15;  # volume in liters
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
  Variable Variable(DNA)
    {
      Value             3300; 
    } 
}


System System( /VacuoleA )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( SIZE )
    {
      Value	 0.03e-15;  # volume in liters
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
      Value     0.7;               
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
}

System System( /VacuoleA/Membrane )
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
}


System System( /VacuoleB )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( SIZE )
    {
      Value	 0.06e-15;  # volume in liters
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
      Value     -0.6;               
    } 
  Variable Variable( ORIGINY )
    {
      Value     0.2;          
    } 
  Variable Variable( ORIGINZ )
    {
      Value     0.1;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
}

System System( /VacuoleB/Membrane )
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
}

System System( /VacuoleC )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( SIZE )
    {
      Value	 0.06e-15;  # volume in liters
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
      Value     -0.5;               
    } 
  Variable Variable( ORIGINY )
    {
      Value     0.5;          
    } 
  Variable Variable( ORIGINZ )
    {
      Value     0.5;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
}

System System( /VacuoleC/Membrane )
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
}

System System( /LysosomeA )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( SIZE )
    {
      Value	 0.06e-15;  # volume in liters
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
      Value     0.55;               
    } 
  Variable Variable( ORIGINY )
    {
      Value     0;          
    } 
  Variable Variable( ORIGINZ )
    {
      Value     -0.55;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
}

System System( /LysosomeA/Membrane )
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
}



System System( /LysosomeB )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( SIZE )
    {
      Value	 0.06e-15;  # volume in liters
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
      Value     0.3;               
    } 
  Variable Variable( ORIGINY )
    {
      Value     -0.6;          
    } 
  Variable Variable( ORIGINZ )
    {
      Value     0;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
}

System System( /LysosomeB/Membrane )
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
}

System System( /LysosomeC )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, SURFACEX, SURFACEY, SURFACEZ)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( SIZE )
    {
      Value	 0.06e-15;  # volume in liters
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
      Value     -0.3;               
    } 
  Variable Variable( ORIGINY )
    {
      Value     -0.55;          
    } 
  Variable Variable( ORIGINZ )
    {
      Value     0;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
}

System System( /LysosomeC/Membrane )
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
}


