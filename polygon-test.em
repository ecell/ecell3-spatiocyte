Stepper SpatiocyteStepper(SS01)
{
    # VoxelRadius 0.005;
    VoxelRadius 0.005;
}

System System(/)
{
    StepperID SS01;

    Variable Variable(TYPE) { Value 0; } 
    Variable Variable(SHAPE) { Value 3; }
	# Variable Variable(VOLUME) {  Value 1.0; }
    Variable Variable(LENGTHX) { Value 1.0; }
    Variable Variable(LENGTHY) { Value 1.0; }
    Variable Variable(LENGTHZ) { Value 1.0; }
    Variable Variable(XYPLANE) { Value 5; }
    Variable Variable(YZPLANE) { Value 5; }
    Variable Variable(XZPLANE) { Value 5; }
    Variable Variable(ORIGINX) { Value 0; } # from -1.0 to +1.0
    Variable Variable(ORIGINY) { Value 0; }
    Variable Variable(ORIGINZ) { Value 0; }
    Variable Variable(VACANT) { Value 0; }

    Variable Variable(A) { Value 0; }
    Variable Variable(B) { Value 0; }
    Variable Variable(C) { Value 0; }

    Process DiffusionProcess(A)
    {
        D 0;
        VariableReferenceList [_ Variable:/:A];
    }

    Process DiffusionProcess(B)
    {
        D 0;
        VariableReferenceList [_ Variable:/:B];
    }

    Process DiffusionProcess(C)
    {
        D 0;
        VariableReferenceList [_ Variable:/:C];
    }

    Process VisualizationLogProcess(loggerMean)
    {
        VariableReferenceList [_ Variable:/:A] [_ Variable:/:B] [_ Variable:/:C];
        LogInterval 1;
    }
  
    Process MoleculePopulateProcess(populate)
    {
        VariableReferenceList [_ Variable:/:A] [_ Variable:/:B] [_ Variable:/:C];
    }
}
