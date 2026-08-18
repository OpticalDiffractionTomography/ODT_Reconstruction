original sources are 
ftp://ftp.wiley.com/public/sci_tech_med/phase_unwrapping/

I modified only 4 files(maskflat.c, util.c, util.h, and maingold.c(the name of which was changed to unwrap2.c)
1. add Free funcion in maskfat.c
2. add void Free(void*); in util.h
3. add MATLAB_MEX_FILE in util.c
4. some modification in maingold.c(unwrap2.c) for the interface to matlab.

