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

#include "DCPDKesslerHigh.h"

#include <cmath>

using namespace std;

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<DCPDKesslerHigh>()
{
  InputParameters params = validParams<Kernel>();
params.addRequiredCoupledVar("v", "Coupled variable");
  return params;
}

DCPDKesslerHigh::DCPDKesslerHigh(const InputParameters & parameters) :
  // You must call the constructor of the base class first
  Kernel(parameters),



_v_var(coupledValue("v")),// _v_var WILL HOLD TEMPERATURE VALUE

    _T_id(coupled("v")) //this is the coupled variable TEMPERATURE identifier which will be used in offDiagQPJacobian

{}

Real DCPDKesslerHigh::computeQpResidual()
{
  // velocity * _grad_u[_qp] is actually doing a dot product
  /* To be copied to QpJacobian */
//Constants
    Real A = pow(10.0,6.649);
Real E = 58160.0;
Real n = 2.192;
Real Kcat = 0.370;
Real R = 8.314;

//Temperature dependence: K(T)
Real KT = A*exp(-E/(R*(_v_var[_qp]+273.0)));

//Cure dependence: f(alpha)
Real falpha = pow(1.0-_u[_qp],n)*(1+Kcat*_u[_qp]);
    /* To be copied to QpJacobian */

//Cure function is product of K(T) and falpha
Real Pc = KT*falpha;


  //cout<<"Pc value is "<<Pc<<endl;
  return  -_test[_i][_qp]*Pc; // changed this, note - sign
//return -_test[_i][_qp]*400.0;
}

Real DCPDKesslerHigh::computeQpJacobian()
{
  // the partial derivative of _grad_u is just _grad_phi[_j]
  /* To be copied to QpJacobian */
//Constants
    Real A = pow(10.0,6.649);
Real E = 58160.0;
Real n = 2.192;
Real Kcat = 0.370;
Real R = 8.314;

//Temperature dependence: K(T)
Real KT = A*exp(-E/(R*(_v_var[_qp]+273.0)));

//Cure dependence: f(alpha)
Real falpha = pow(1.0-_u[_qp],n)*(1+Kcat*_u[_qp]);
    /* To be copied to QpJacobian */

Real dfalpha = pow(1.0-_u[_qp],n)*Kcat - n*pow(1-_u[_qp],n-1)*(1+Kcat*_u[_qp]);
//Cure function is product of K(T) and falpha
Real Pc = KT*falpha;

  Real dPc = KT*dfalpha;
    //cout<<"dPcbydalpha value is "<<dPc<<endl;

  return -_test[_i][_qp]*dPc*_phi[_j][_qp];
  //return 0.0;
  //note that chain rule is used not in dPc calculation, but used in return statement by multiplying with _phi[_j][_qp]
  //see MOOSE manual on how to calculate Jacobian
  //I changed this, and note the - sign

}

//off diagonal Jacobian needs to be added !!!!!!!!!!
//
Real DCPDKesslerHigh::computeQpOffDiagJacobian(unsigned int jvar)
{
    cout<<"offDiagonalJacobian executed"<<endl;
    if (jvar==_T_id)
    {
  /* To be copied to QpJacobian */
//Constants
    Real A = pow(10.0,6.649);
Real E = 58160.0;
Real n = 2.192;
Real Kcat = 0.370;
Real R = 8.314;

//Temperature dependence: K(T)
Real KT = A*exp(-E/(R*(_v_var[_qp]+273.0)));

//Cure dependence: f(alpha)
Real falpha = pow(1.0-_u[_qp],n)*(1+Kcat*_u[_qp]);
    /* To be copied to QpJacobian */

Real dKT = KT*(E/R)*1/pow(_u[_qp],2);

Real dPcbydT=dKT*falpha;

    cout<<"dPcbydT value is "<<dPcbydT<<endl;

    return -_test[_i][_qp]*dPcbydT*_phi[_j][_qp];
     //return 0.0;
     //_v_var[_qp] from residual becomes _phi[_j][_qp] in the off diagonal jacobian
    }

    return 0.0;
}
