#FC in the filename means Fully Coupled

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 400  #400 elements is the base number
  ny = 1
  xmax = 1 # non dimensional length of test chamber
  ymax = .01 # Test chamber radius
  
[]

[Variables]
  active = 'Temperature Cure'

  [./Temperature]
    order = FIRST
    family = LAGRANGE
    initial_condition = '0'
  [../]
  
  [./Cure]
    order = FIRST
    family = LAGRANGE
    initial_condition = '0.0000'
  [../]
[]

[Kernels]
  active = 'tempdiff tempderv coupledcurederv cureformula curederv' 

  [./tempdiff] #temperature diffusion term
    type = TempDiffusion
    variable = Temperature
    
  [../]

  [./coupledcurederv]
    type = CoupledCureTimeDerivative
    variable = Temperature
    v = Cure
    
  [../]
  
  [./tempderv]
  type = HeatConductionTimeDerivative
  variable = Temperature
  lumping='false'
  
  [../]
  
  [./cureformula]
  type=DCPDnonD
  variable = Cure
  v = Temperature #this is the coupled variable
  Ttrig = '300'
  Tintl = '20'
  [../]
  
  [./curederv]
  type=TimeDerivative
  variable = Cure
  lumping='false'
  [../]  
  
[]

[BCs]
  active = 'templeft'

  [./templeft]
    type = DirichletBC
    #type = NeumannBC
    variable = Temperature
    boundary = 'left'
    value = '1'
  [../]
  
  [./tempright]
    type = DirichletBC
    #type = NeumannBC
    variable = Temperature
    boundary = 'right'
    value = '30'
  [../]
          
[]

[Materials]
  [./polymer]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'specific_heat density Hr TConductivity'
    prop_values = '1.0 1.0 1.1409 3.5245e-10' # polymer attributes
  [../]
[]

#[Preconditioning]  #to be removed in an ordinary run
#  [./SMP]
#    type = FDP
#    full = true
#  [../]
#[]

[Executioner]
  type = Transient
  num_steps = 10 ####################################
  
  [./TimeStepper]
    type = ConstantDT    
    dt = 5729.6
  [../]
  
  [./TimeIntegrator]
   type = ImplicitEuler
[../]
  
  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

[]


[Outputs] 
  execute_on = 'timestep_end' # Limit the output to timestep end (removes initial condition)
  [./console]
    type = Console
    perf_log = true          # enable performance logging
  [../]
  [./exodus]
    type = Exodus
    execute_on = 'initial timestep_end' # output the initial condition for the file
    file_base = adaptoff  # set the file base (the extension is automatically applied), contains only the refine and coarsen parameter values 
    interval = 10           # only output every 10 step
  [../]  
[]


