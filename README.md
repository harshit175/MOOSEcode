

# MOOSEcode

The code solves the coupled thermochemical PDEs for the  frontal polymerisation problem . The equations are solved in the non dimensional form. The project was done while pursuing MEng in Mechanical Engineering at UIUC.  Currently the documentation is not in a form readable by general public. 

# Short Description
The objective of the project was to model the frontal polymerization process.  It was being done in collaboration with a group in the chemistry department working on experiments.  The setup is a microchannel of around 1 mm diameter containing a monomer.  The microchannel is surrounded by an inert matrix for insulation.  The polymerization progresses in the microchannel as a self sustaining thermal front upon ignition with a soldering iron from one end.  Thermochemical equations governing the process were numerically solved by developing a finite element code.  Parameters like front speed, acceleration and temperature profile were computed.  The application of this phenomenon is in rapid manufacturing of composites.  The coding was done an open source finite element framework known as MOOSE.

# Code Details for the kernels
The TempDiffusion and CoupledCureTimeDerivative kernels were originally coded for the  dimensional  form,  but they are used for the non dimensional form of the equation, this is done by calculating the non-dimensional parameters and entering them lumped, under "Tconductivity" and "Hr". The density and specific heat in are each entered as 1.  The time step used is the non dimensionalized timestep. (Tao=A.t)

The src folder contains different kernels for the cure kinetics, like CureformulaS1, Cure FormulaS2 etc. the nondimensionalized formuala  is coded in DCPDnonD which  is  used in the input file adapton.i and adaptoff.i.  All other cure formulas are  from earlier work.

The Kernels&Matprop.pdf file shows the dimensional equations and the non-dimensionalized equations with the kernels labelled. The DimensionalandNonDimensional.pdf goes into greater detail of the formulation.
