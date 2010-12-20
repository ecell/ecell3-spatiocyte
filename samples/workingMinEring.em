# A test model for Spatiocyte 
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 3e-8;
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
      Value     1;      # { 0: Spherical (uses SIZE) 
                        #   1: Rod (uses SIZE, LENGTHY == 2*radius)
                        #   2: Cubic (uses SIZE, PLANEYZ, PLANEXZ, PLANEXY)
                        #   3: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                        #   4: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable( LENGTHX )
    {
      Value     0e-6;        # in meters
    } 
  Variable Variable( LENGTHY )
    {
      Value     1e-6;        # in meters
    } 
  Variable Variable( LENGTHZ )
    {
      Value     0e-6;        # in meters
    } 

  Variable Variable( SIZE )
    {
      Value	 2.88e-15;   # volume in liters
                              # RAM (memory) required volume/
                              # 0.85e-14 liter max for 8GB memory
                              # e.coli sphere radius: 0.5e-6 m 
                              # e.coli cylinder length: 3e-6 m
      # e.coli volume = (4*PI*(0.5e-6)^3)/3 + PI*(0.5e-6)^(2)*(3e-6)
      #               = 2.88e-18 m^3
      #               = 2.88e-15 l
    } 
  Variable Variable( PLANEYZ )
    {
      Value     0;               # { 0: Reflective (Impermeable) 
                                 #   1: Periodic (Permeable) }
    } 
  Variable Variable( PLANEXZ )
    {
      Value     0;               # { 0: Reflective (Impermeable) 
                                 #   1: Periodic (Permeable) }
    } 
  Variable Variable( PLANEXY )
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
  Variable Variable( MinD_ATP )
    {
      Value             500; 
    } 
  Variable Variable( MinD_ADP )
    {
      Value             500; 
    } 
  Variable Variable( MinE )
    {
      Value             700; 
    } 
  Process DiffusionProcess( diffuse )
    {
      VariableReferenceList [ _ Variable:/:MinD_ATP ]
                            [ _ Variable:/:MinD_ADP ]
                            [ _ Variable:/:MinE ];
      D 2.5e-12;
    }
  Process VisualizationLogProcess( loggerMean )
    {
      VariableReferenceList [ _ Variable:/Surface:MinD_ATP_m_MinE ]
                            [ _ Variable:/Surface:MinD_ADP_m_MinE ]
                            [ _ Variable:/Surface:MinE_m ]
                            [ _ Variable:/Surface:MinD_ATP_m ];
      LogInterval 0.05;
      FileName "visualLog0.dat";
    }
  # MinD_ADP -> MinD_ATP
  Process SpatiocyteNextReactionProcess( MinD_ADPtoMinD_ATP )
    {
      VariableReferenceList [_ Variable:/:MinD_ADP -1]
                            [_ Variable:/:MinD_ATP 1];
      k 6; # in s^{-1}
    }
  Process MoleculePopulateProcess( populate )
    {
      VariableReferenceList [_ Variable:/:MinD_ATP ]
                            [_ Variable:/:MinD_ADP ]
                            [_ Variable:/:MinE ]
                            [_ Variable:/Surface:MinD_ATP_m ]
                            [_ Variable:/Surface:MinD_ATP_m_MinE ]
                            [_ Variable:/Surface:MinD_ADP_m_MinE ]
                            [_ Variable:/Surface:MinE_m ];
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

  Variable Variable( MinD_ATP_m )
    {
      Value             0; 
    } 
  Variable Variable( MinD_ATP_m_MinE )
    {
      Value             0; 
    } 
  Variable Variable( MinD_ADP_m_MinE )
    {
      Value             0; 
    } 
  Variable Variable( MinE_m )
    {
      Value             0; 
    } 
  Process DiffusionProcess( diffuseS )
    {
      VariableReferenceList [ _ Variable:/Surface:MinD_ATP_m ];
      D 1e-14;
    }
  Process DiffusionProcess( diffusee )
    {
      VariableReferenceList [ _ Variable:/Surface:MinE_m ];
      D 1e-14;
    }
  Process DiffusionProcess( diffuseSss )
    {
      VariableReferenceList [ _ Variable:/Surface:MinD_ATP_m_MinE ];
      D 1e-14;
    }
  Process DiffusionProcess( diffuseSs )
    {
      VariableReferenceList [ _ Variable:/Surface:MinD_ADP_m_MinE ];
      D 1e-14;
    }
  # MinD_ATP + LIPID -> MinD_ATP_m
  Process DiffusionInfluencedReactionProcess( MinD_ATP_2_MinD_ATP_m ) 
    {
      VariableReferenceList   [ _ Variable:/Surface:LIPID -1 ]
                              [ _ Variable:/:MinD_ATP -1 ]
                              [ _ Variable:/Surface:MinD_ATP_m 1 ];
      k 0.0125e-6; # in ms^{-1}
      # In 2006.fange.plos.comput.biol the value of k is
      # 0.0125 ums^{1}
    } 
  # MinD_ATP_m + MinD_ATP -> MinD_ATP_m + MinD_ATP_m
  Process DiffusionInfluencedReactionProcess( MinD_ATP_mtoMinD_ATP_m_MinE )
    {
      VariableReferenceList   [ _ Variable:/:MinD_ATP -1 ]
                              [ _ Variable:/Surface:MinD_ATP_m -1 ]
                              [ _ Variable:/Surface:MinD_ATP_m 1 ]
                              [ _ Variable:/Surface:MinD_ATP_m 1 ];
      k 1.49448498e-20; # in m^3s^{-1}
      # In 2006.fange.plos.comput.biol the value of k is
      # 9e+6 M^{-1}s^{-1} = 9e+6*1e-3/6.0221415e+23 m^3s^{-1}
      #                   = 1.49448498e-20 m^3s^{-1}
    } 
  # MinD_ATP_m_MinE -> MinD_ADP + MinE
  Process DiffusionInfluencedReactionProcess( MinD_ATP_m_MinEtoMinD_ATP_m_MinE )
    {
      VariableReferenceList   [ _ Variable:/Surface:MinD_ATP_m -1 ]
                              [ _ Variable:/:MinE -1 ]
                              [ _ Variable:/Surface:MinD_ATP_m_MinE 1 ];
      k 9.23259608e-20;
      #p 1;
      # k 9.23259608e-20; # in m^3s^{-1}
      # In 2006.fange.plos.comput.biol the value of k is
      # 5.56e+7 M^{-1}s^{-1} = 5.6e+7*1e-3/6.0221415e+23 m^3s^{-1}
      #                      = 9.223259608-20 m^3s^{-1}
    } 
  Process SpatiocyteNextReactionProcess( MinD_ATP_m_MinEtoMinD_ADP_MinE )
    {
      VariableReferenceList   [ _ Variable:/Surface:MinD_ATP_m_MinE -1 ]
                              [ _ Variable:/:MinD_ADP 1 ]
                              [ _ Variable:/:MinE 1 ];
      k 0.5; # in s^{-1}
    } 
  Process DiffusionInfluencedReactionProcess( polymerize2s )
    {
      VariableReferenceList   [ _ Variable:/Surface:MinD_ATP_m_MinE -1 ]
                              [ _ Variable:/Surface:MinD_ATP_m_MinE -1 ]
                              [ _ Variable:/Surface:MinD_ATP_m_MinE 1 ]
                              [ _ Variable:/Surface:MinD_ADP_m_MinE 1 ];
      k 9.23259608e-15;
    } 
  Process SpatiocyteNextReactionProcess( MinD_ATP_mp22s )
    {
      VariableReferenceList   [ _ Variable:/Surface:MinD_ADP_m_MinE -1 ]
                              [ _ Variable:/Surface:MinE_m 1 ]
                              [ _ Variable:/:MinD_ADP 1 ];
      k 0.1; # in m^3s^{-1}
    } 
  Process SpatiocyteNextReactionProcess( MinD_ATP_mp2s )
    {
      VariableReferenceList   [ _ Variable:/Surface:MinE_m -1 ]
                              [ _ Variable:/:MinE 1 ];
      k 0.1; # in m^3s^{-1}
    } 
}




