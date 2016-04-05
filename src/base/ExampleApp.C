/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#include "ExampleApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

// Example 2 Includes
//#include "ExampleConvection.h"           // <- New include for our custom kernel
#include "TempDiffusion.h"               // <- I included TempDiffusion kernel
#include "HeatConductionTimeDerivative.h"// <- I included time derivative kernel
#include "CureFormula.h" // I included
#include "CureFormulaDCPD.h" // I included for DCPD medium conc. from Kessler
#include "DCPDKesslerLow.h" // I included for DCPD low conc. from Kessler
#include "DCPDKesslerHigh.h" // I included for DCPD high conc. from Kessler
#include "CureFormulaACR.h" // I included for Acrylate
#include "CureFormulaS1.h" // I included
#include "CureFormulaS2.h" // I included
#include "CureFormulaS3.h" // I included
#include "DCPDnonD.h" //I included for non dimensional, DCPD med conc. from Kessler
#include "DCPDnonDgeneral.h" //I included for non dimensional, DCPD general conc. from Kessler
#include "CoupledCureTimeDerivative.h"
//#include "TimeDerivative.h"
//#include "CoupledTimeDerivative.h"
template<>
InputParameters validParams<ExampleApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

ExampleApp::ExampleApp(InputParameters parameters) :
    MooseApp(parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ExampleApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ExampleApp::associateSyntax(_syntax, _action_factory);
}

ExampleApp::~ExampleApp()
{
}

void
ExampleApp::registerObjects(Factory & factory)
{
  // Register any custom objects you have built on the MOOSE Framework
  registerKernel(CoupledCureTimeDerivative);  // <- registration
  registerKernel(TempDiffusion);      // <- registration for temperature diffusion created by me
  registerKernel(HeatConductionTimeDerivative);     // <- registration for heat conduction time derivative kernel
  //registerKernel(CoupledTimeDerivative); //will be used for cure time derivative
  //registerKernel(TimeDerivative);
   registerKernel(CureFormula); //contains the Frulloni cure function
   registerKernel(CureFormulaDCPD); //contains the DCPD medium conc. cure function from Kessler paper
    registerKernel(DCPDKesslerLow); //contains the DCPD low conc. cure function from Kessler paper
   registerKernel(DCPDKesslerHigh); //contains the DCPD high conc. cure function from Kessler paper
    registerKernel(DCPDnonD); //contains DCPD medium conc. nondimensional cure function from Kessler paper
registerKernel(DCPDnonDgeneral); //contains DCPD general nondimensional cure function from Kessler paper
registerKernel(CureFormulaACR);  //contains the cure formula I arrived at from Ian's EP-21 data
registerKernel(CureFormulaS1); //contains the simple cure function
registerKernel(CureFormulaS2); //contains the simple cure function with linear temperature dependence
registerKernel(CureFormulaS3); //contains the simple cure function with arrhenius temperature dependence

}

void
ExampleApp::registerApps()
{
  registerApp(ExampleApp);
}

void
ExampleApp::associateSyntax(Syntax& /*syntax*/, ActionFactory & /*action_factory*/)
{
}
