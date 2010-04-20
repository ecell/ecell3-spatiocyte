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
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.5e-6;        # in meters
    } 
  Variable Variable( SIZE )
    {
      Value	 5e-13;  # volume in liters
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
  Process VisualizationLogProcess(loggerMean)
    {
      VariableReferenceList [_ Variable:/CellB/Membrane:LIPID]
                            [_ Variable:/CellA:A]
                            [_ Variable:/CellA:B]
                            [_ Variable:/CellA:C]
                            [_ Variable:/CellA/Membrane:D]
                            [_ Variable:/CellA/Membrane:E]
                            [_ Variable:/CellA/Membrane:F]
                            [_ Variable:/CellA/Membrane:LIPID]
                            [_ Variable:/CellA/Nucleus:G]
                            [_ Variable:/CellA/Nucleus:H]
                            [_ Variable:/CellA/Nucleus:I]
                            [_ Variable:/CellA/Nucleus/Membrane:J]
                            [_ Variable:/CellA/Nucleus/Membrane:K]
                            [_ Variable:/CellA/Nucleus/Membrane:L]
                            [_ Variable:/CellA/Nucleus/Nucleolus:DNA]
                            [_ Variable:/CellA/Nucleus/Membrane:LIPID]
                            [_ Variable:/CellA/MitochondrionA/Membrane:LIPID]
                            [_ Variable:/CellA/MitochondrionB/Membrane:LIPID]
                            [_ Variable:/CellA/MitochondrionC/Membrane:LIPID]
                            [_ Variable:/CellA/VacuoleA/Membrane:LIPID]
                            [_ Variable:/CellA/VacuoleB/Membrane:LIPID]
                            [_ Variable:/CellA/VacuoleC/Membrane:LIPID]
                            [_ Variable:/CellA/LysosomeA/Membrane:LIPID]
                            [_ Variable:/CellA/LysosomeB/Membrane:LIPID]
                            [_ Variable:/CellA/LysosomeC/Membrane:LIPID];
      LogInterval 0.05;
      FileName "visualLog0.dat";
    }
}

System System( /CellB )
{
  StepperID       SS; 
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.5e-6;        # in meters
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
      Value     0.46;       
    } 
  Variable Variable( VACANT )
    {
      Value             0; 
    } 
}

System System( /CellB/Membrane )
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

System System( /CellA )
{
  StepperID       SS; 
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.5e-6;        # in meters
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
      Value     -0.46;       
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
}

System System( /CellA/Membrane )
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

System System( /CellA/Nucleus )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.5e-6;        # in meters
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

System System( /CellA/MitochondrionA )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     1;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.2e-6;        # in meters
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

System System( /CellA/MitochondrionA/Membrane )
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

System System( /CellA/MitochondrionB )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     1;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.2e-6;        # in meters
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

System System( /CellA/MitochondrionB/Membrane )
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

System System( /CellA/MitochondrionC )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     1;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.2e-6;        # in meters
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

System System( /CellA/MitochondrionC/Membrane )
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

System System( /CellA/Nucleus/Membrane )
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


System System( /CellA/Nucleus/Nucleolus )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.2e-6;        # in meters
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


System System( /CellA/VacuoleA )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.2e-6;        # in meters
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

System System( /CellA/VacuoleA/Membrane )
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


System System( /CellA/VacuoleB )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.2e-6;        # in meters
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

System System( /CellA/VacuoleB/Membrane )
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

System System( /CellA/VacuoleC )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.2e-6;        # in meters
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

System System( /CellA/VacuoleC/Membrane )
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

System System( /CellA/LysosomeA )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.2e-6;        # in meters
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

System System( /CellA/LysosomeA/Membrane )
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



System System( /CellA/LysosomeB )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.2e-6;        # in meters
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

System System( /CellA/LysosomeB/Membrane )
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

System System( /CellA/LysosomeC )
{
  StepperID SS;
  Variable Variable( TYPE )
    {
      Value     0;              # { 0: Volume (requires SHAPE, RADIUS, SIZE)
                                #   1: Surface }
    } 
  Variable Variable( SHAPE )
    {
      Value     0;              # { 0: Spherical
                                #   1: Rod (uses RADIUS)
                                #   2: Cubic (uses SURFACEX,SURFACEY,SURFACEZ) }
    } 
  Variable Variable( RADIUS )
    {
      Value     0.2e-6;        # in meters
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

System System( /CellA/LysosomeC/Membrane )
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

