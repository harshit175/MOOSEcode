#FC in the filename means Fully Coupled
[Mesh]
  file=2domain6mmMATLAB.inp    
  uniform_refine = 3
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
    initial_condition = '0.000'
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
  #specific_heat = '1'
  #density = '1'
  #the specific heat and density are added in materials and not in this block
  [../]
  
  [./cureformula]
  type=DCPDnonDgeneral
  variable = Cure
  v = Temperature #this is the coupled variable
  Ttrig = '339.44'
  Tintl = '20'
  _E = '51100'
  _n = '1.9270'
  _Kcat = '0.3650'
  
  
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
    boundary = 'a_chanbot'
    value = '0.85'
    
  [../]
        
#[]

[Problem]  #makes the problem axisymmetric
  type = FEProblem
  coord_type = RZ
  rz_coord_axis = Y
[]

#length is 5.5*33mm, Rm is 5.5mm
[Materials]
  [./DCPD]
    type = GenericConstantMaterial
    block = '4'
    prop_names = 'specific_heat density Hr TConductivity'
    prop_values = '1.0 1.0 1.000 1.8644e-8' # polymer attributes
  [../]
  
  [./PDMS]
    type = GenericConstantMaterial
    block = '5'
    prop_names = 'specific_heat density Hr TConductivity'
    prop_values = '1 .8029 0 3.3117e-8' # PDMS attributes
  [../]
  
  #[./Copper]
  #  type = GenericConstantMaterial
  #  block = '4'
  #  prop_names = 'specific_heat density Hr TConductivity'
  #  prop_values = '1 2.4399 0 7.3474e-5' # wire attributes
  #[../]

[]

#[Preconditioning]  #to be removed in an ordinary run
#  [./SMP]
#    type = FDP
#    full = true
#  [../]
#[]

[Executioner]
  type = Transient
  num_steps = 12000 
  #l_rel_tol = 1e-7
  nl_rel_tol = 1e-6
  
  [./TimeStepper]
    type = ConstantDT    
    dt = 11459.2 
  [../]
  
  [./TimeIntegrator]
   type = ImplicitEuler
[../]
  
  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  
[]

[Adaptivity]
  marker = errorfrac  #this line when commented, switches off adaptivity
  max_h_level = 3
  #steps = 2 #this line gets ignored in a transient run
  [./Indicators]
    [./error]
      type = GradientJumpIndicator
      variable = Temperature
      outputs = none
    [../]
    
  [../]
  [./Markers]
    [./errorfrac]
      type = ErrorFractionMarker
      refine = 0.9
      coarsen = 0.02
      indicator = error
      outputs = none
    [../]
    
  [../]
[]


[Controls]
  [./bcs]
    type = TimePeriod
    disable_objects = 'BCs::templeft'
    start_time = '1718880' #'2005360'  '1317808'   
    execute_on = 'initial timestep_begin'
  [../]  
[]

[Outputs] 
  execute_on = 'timestep_end' # Limit the output to timestep end (removes initial condition)
  [./console]
    type = Console
    perf_log = true          # enable performance logging    
  [../]
  [./cp]
    type = Checkpoint
    num_files = 4
    interval = 1
  [../]
  [./exodus]
    type = Exodus
    execute_on = 'initial timestep_end' # output the initial condition for the file
    file_base = Rc5_5mm # set the file base (the extension is automatically applied) 
    interval = 10           # only output every 10 step
  [../]  
[]

