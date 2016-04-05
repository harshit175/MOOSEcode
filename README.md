# MOOSEcode

The code solves the coupled thermochemical PDEs for the  frontal polymerisation problem . The equations are solved in the non dimensional form.

The TempDiffusion and CoupledCureTimeDerivative kernels were originally coded for the  dimensional  form,  but they are used for the non dimensional form of the equation, this is done by calculating the non-dimensional parameters and entering them lumped, under "Tconductivity" and "Hr". The density and specific heat in are each entered as 1.  The time step used is the non dimensionalized timestep. (Tao=A.t)

The src folder contains different kernels for the cure kinetics, like CureformulaS1, Cure FormulaS2 etc. the nondimensionalized formuala  is coded in DCPDnonD which  is  used in the input file adapton.i and adaptoff.i.  All other cure formulas are  from earlier work.
