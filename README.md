#DEBUGINFO

The code runs fine without adaptivity (using adaptoff.i) in the debug mode but gives a segfault when run the usual way.

The code runs till first timestep in the debug mode, giving an arithmetic exception (SIGFPE), but gives a segfault in the debug mode without progressing at all on the time stepping.

The above observations are strange.  

# MOOSEcode

The code solves the coupled thermochemical PDEs for the  frontal polymerisation problem . The equations are solved in the non dimensional form. The project forms a part of my PhD thesis 

The TempDiffusion and CoupledCureTimeDerivative kernels were originally coded for the  dimensional  form,  but they are used for the non dimensional form of the equation, this is done by calculating the non-dimensional parameters and entering them lumped, under "Tconductivity" and "Hr". The density and specific heat in are each entered as 1.  The time step used is the non dimensionalized timestep. (Tao=A.t)

The src folder contains different kernels for the cure kinetics, like CureformulaS1, Cure FormulaS2 etc. the nondimensionalized formuala  is coded in DCPDnonD which  is  used in the input file adapton.i and adaptoff.i.  All other cure formulas are  from earlier work.

The Kernels&Matprop.pdf file shows the dimensional equations and the non-dimensionalized equations with the kernels labelled. The DimsnsionalandNonDimensional.pdf goes into greater detail of the formulation.
