*************************************************************************
README
2017.09.11
*************************************************************************
the first argument stands for the radius[um] of the bead.
the second argument stands for the concentration[g/g].
e.g.
     ./mie 1 0.025
     send 1 as the radius of the bead,0.025 as concentration and output:

     "input_data/input.txt"           mua,mus
     "input_data/index.txt"           index of water, index of bead
     "phase/phaseOutput_xxxnm.txt"    phase function of each wavelength

     these output file is place carefully for MCML_Phantom program to read.
     if changing the path of these file, check the input file path of
     MCML_Phantom program either.

If changing other parameters is needed,you should modify the MisStruct object.

modified by SHI-CHENG, TU 2017.09.11
