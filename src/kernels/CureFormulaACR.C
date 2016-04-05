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

#include "CureFormulaACR.h"

#include <cmath>

using namespace std;

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<CureFormulaACR>()
{
  InputParameters params = validParams<Kernel>();
params.addRequiredCoupledVar("v", "Coupled variable");
  return params;
}

CureFormulaACR::CureFormulaACR(const InputParameters & parameters) :
  // You must call the constructor of the base class first
  Kernel(parameters),



_v_var(coupledValue("v")),// _v_var WILL HOLD TEMPERATURE VALUE

    _T_id(coupled("v")) //this is the coupled variable TEMPERATURE identifier which will be used in offDiagQPJacobian

{}

Real CureFormulaACR::computeQpResidual()
{
 
  /* To be copied to QpJacobian */
Real a=1.0;
 

  Real k02=a*exp(44.2);
  
  Real Eatt2=127038.0;
  Real mm=0.1127;
  Real nn=1.6563;
  Real R=8.314;
  

  Real k2=k02*exp(-Eatt2/(R*(_v_var[_qp]+273.0)));

  Real p = .02184;
  Real q = -.3107;
  Real alphamax = 1.0;
  if ( _v_var[_qp]<60.0)
	   alphamax = p*_v_var[_qp]+q;
    /* To be copied to QpJacobian */


  Real Pc=(k2*pow(_u[_qp],mm))*pow((alphamax-_u[_qp]),nn);
  //cout<<"Pc value is "<<Pc<<endl;
  return  -_test[_i][_qp]*Pc; // changed this, note - sign
//return -_test[_i][_qp]*400.0;
}

Real CureFormulaACR::computeQpJacobian()
{
  // the partial derivative of _grad_u is just _grad_phi[_j]
/* To be copied to QpJacobian */
Real a=1.0;
 

  Real k02=a*exp(44.2);
  
  Real Eatt2=127038.0;
  Real mm=0.1127;
  Real nn=1.6563;
  Real R=8.314;
  

  Real k2=k02*exp(-Eatt2/(R*(_v_var[_qp]+273.0)));

  Real p = .02184;
  Real q = -.3107;
  Real alphamax = 1.0;
  if (_v_var[_qp]<60.0)
	   alphamax = p*_v_var[_qp]+q;
    /* To be copied to QpJacobian */

  Real Pc=(k2*pow(_u[_qp],mm))*pow((alphamax-_u[_qp]),nn);
  Real dPc = (k2*mm*pow(_u[_qp],mm-1.0)*pow((alphamax-_u[_qp]),nn)-(k2*pow(_u[_qp],mm))*nn*pow((alphamax-_u[_qp]),nn-1.0));
    //cout<<"dPcbydalpha value is "<<dPc<<endl;

  return -_test[_i][_qp]*dPc*_phi[_j][_qp];
  //return 0.0;
  //note that chain rule is used not in dPc calculation, but used in return statement by multiplying with _phi[_j][_qp]
  //I changed this, and note the - sign

}

//off diagonal Jacobian needs to be added !!!!!!!!!!
//
Real CureFormulaACR::computeQpOffDiagJacobian(unsigned int jvar)
{
    cout<<"offDiagonalJacobian executed"<<endl;
    if (jvar==_T_id)
    {
/* To be copied to QpJacobian */
Real a=1.0;
 

  Real k02=a*exp(44.2);
  
  Real Eatt2=127038.0;
  Real mm=0.1127;
  Real nn=1.6563;
  Real R=8.314;
  

  Real k2=k02*exp(-Eatt2/(R*(_v_var[_qp]+273.0)));

  Real p = .02184;
  Real q = -.3107;
  Real alphamax = 1.0;
  if (_v_var[_qp]<60.0)
	   alphamax = p*_v_var[_qp]+q;
    /* To be copied to QpJacobian */


   
    Real dk2bydT=k2*(Eatt2/R)*pow((1.0/(_v_var[_qp]+273.0)),2.0);

    Real dPcbydT=(dk2bydT*pow(_u[_qp],mm))*pow((alphamax-_u[_qp]),nn);
    cout<<"dPcbydT value is "<<dPcbydT<<endl;

    return -_test[_i][_qp]*dPcbydT*_phi[_j][_qp];
     //return 0.0;
     //_v_var[_qp] from residual becomes _phi[_j][_qp] in the off diagonal jacobian
    //note that chain rule is applied in the return statement
    }

    return 0.0;
}
