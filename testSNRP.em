# A test model for Spatiocyte 
# written by Satya Arjunan <satya.arjunan(a)gmail.com>

Stepper SpatiocyteStepper(SS)
{
  VoxelRadius 4.4e-9;
}

Stepper ODEStepper(DE)
{
}

System System( / )
{
  StepperID       SS; 
  Variable Variable(GEOMETRY)
    {
      Value 3;         # { 0: Cuboid (uses LENGTHX, LENGTHY, LENGTHZ)
                       #   1: Ellipsoid (uses LENGTHX, LENGTHY, LENGTHZ) }
                       #   2: Cylinder (uses LENGTHX, LENGTHY=2*radius)
                       #   3: Rod (uses LENGTHX, LENGTHY=2*radius)
                       #   4: Torus (uses LENGTHX, LENGTHY, LENGTHZ) }
                       #   5: Pyramid (uses LENGTHX, LENGTHY, LENGTHZ) }
    } 
  Variable Variable(LENGTHX)
    {
      Value 4.5e-6;      # m
    } 
  Variable Variable(LENGTHY)
    {
      Value 1e-6;      # m
    } 
  Variable Variable(VACANT)
    {
      Value 0; 
    } 
  Variable Variable(MinD1)
    {
      Value 0;         # molecule number
    } 
  Variable Variable(MinD2)
    {
      Value 0;         # molecule number
      Name "HD";
    } 
  Variable Variable(MinD3)
    {
      Value 0;         # molecule number
      Name "HD";
    } 
  Variable Variable(MinD4)
    {
      Value 0;         # molecule number
      Name "HD";
    } 
}

System System(/Surface)
{
  StepperID SS;

  Variable Variable(DIMENSION)
    {
      Value 2;         # { 3: Volume
                       #   2: Surface }
    }
  Variable Variable(VACANT)
    {
      Value 0;
    }

  Variable Variable(MinD1)
    {
      Value 2000;         # molecule number
    }
  Process DiffusionProcess( diff )
    {
      VariableReferenceList [_ Variable:/Surface:MinD1]
                            [_ Variable:/:MinD1];
    }
  Process MassActionProcess( ES_to_E_P )
    {
      StepperID DE;
      VariableReferenceList [_ Variable:/Surface:MinD1 -1]
                            [_ Variable:/:MinD1 1];
      k 2.3e-3;
    }  

  Variable Variable(MinD2)
    {
      Value 2000;         # molecule number
      Name "HD";
    }
  Process SpatiocyteNextReactionProcess(sph)
    {
      VariableReferenceList   [_ Variable:/Surface:MinD2 -1]
                              [_ Variable:/:MinD2 1];
      k 2.3e-3;
    }  

  Variable Variable(MinD3a)
    {
      Value 1000;         # molecule number
    }
  Variable Variable(MinD3b)
    {
      Value 1000;         # molecule number
    }
  Process MoleculePopulateProcess(populate3)
    {
      VariableReferenceList  [_ Variable:/Surface:MinD3a]
                             [_ Variable:/Surface:MinD3b];
    }
  Process SpatiocyteNextReactionProcess(sp)
    {
      VariableReferenceList   [_ Variable:/Surface:MinD3a -1]
                              [_ Variable:/:MinD3 1];
      k 2.3e-3;
    }  
  Process SpatiocyteNextReactionProcess(sp2)
    {
      VariableReferenceList   [_ Variable:/Surface:MinD3b -1]
                              [_ Variable:/:MinD3 1];
      k 2.3e-3;
    }  

  Variable Variable(MinDpe)
    {
      Value 1000;         # molecule number
    }
  Variable Variable(MinDpg)
    {
      Value 1000;         # molecule number
    }
  Variable Variable(PE)
    {
      Value 121000;         # molecule number
    }
  Variable Variable(PG)
    {
      Value 121000;         # molecule number
    }
  Process MoleculePopulateProcess(populate4)
    {
      VariableReferenceList  [_ Variable:/Surface:PE]
                             [_ Variable:/Surface:PG];  
    }
  Process MoleculePopulateProcess(populate5)
    {
      VariableReferenceList  [_ Variable:/Surface:MinDpe]
                             [_ Variable:/Surface:MinDpg];  
    }
  Process DiffusionProcess(diffpe)
    {
      VariableReferenceList   [_ Variable:/Surface:MinDpe]
                              [_ Variable:/Surface:PE -1];
      D 0.02e-12; 
    }
  Process DiffusionProcess(diffpg)
    {
      VariableReferenceList   [_ Variable:/Surface:MinDpg]
                              [_ Variable:/Surface:PG -1];
      D 0.02e-12; 
    }
  Process SpatiocyteNextReactionProcess(pg2)
    {
      VariableReferenceList   [_ Variable:/Surface:MinDpg -1]
                              [_ Variable:/:MinD4 1];
      k 2.3e-3;          # s^{-1}
    }  
  Process SpatiocyteNextReactionProcess(pe2)
    {
      VariableReferenceList   [_ Variable:/Surface:MinDpe -1]
                              [_ Variable:/:MinD4 1];
      k 2.3e-3;          # s^{-1}
    }  
  }
