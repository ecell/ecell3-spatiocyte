Stepper SpatiocyteStepper(SS) { VoxelRadius 6e-9; }
Stepper ODEStepper(DE) { MaxStepInterval 1e-5; }

System System(/)
{
  StepperID       SS;
  Variable Variable(GEOMETRY)
    {
      Value 5;
    }
  Variable Variable(LENGTHX)
    {
      Value 1e-6;
    }
  Variable Variable(LENGTHY)
    {
      Value 1e-6;      
    }
  Variable Variable(LENGTHZ)
    {
      Value 1e-6;     
    }
  Variable Variable(VACANT)
    {
      Value 0;
    }
    Variable Variable( A )
    {
        Name HD;
        Value   10000;
    }
    Variable Variable( B )
    {
        Name HD;
        Value   10000;
    }
    Variable Variable( C )
    {
        Name HD;
        Value  0;  
    }
    Variable Variable( Am )
    {
        Name HD;
        Value   10000;
    }
    Variable Variable( Bm )
    {
        Name HD;
        Value   10000;
    }
    Variable Variable( Cm )
    {
        Name HD;
        Value  0;  
    }
   Process SpatiocyteNextReactionProcess( reaction )
    {
      VariableReferenceList [_ Variable:/:A -2]   
                            [_ Variable:/:B -1]   
                            [_ Variable:/:C 1];    
      k                     1e-45;
    }
   Process MassActionProcess( reaction2 )
    {
      StepperID       DE;
      VariableReferenceList [_ Variable:/:Am -2]   
                            [_ Variable:/:Bm -1]   
                            [_ Variable:/:Cm 1];    
      k                     1e-45;
    }
}
