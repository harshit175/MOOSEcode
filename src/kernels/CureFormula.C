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

#include "CureFormula.h"

#include <cmath>

using namespace std;

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<CureFormula>()
{
  InputParameters params = validParams<Kernel>();
params.addRequiredCoupledVar("v", "Coupled variable");
  return params;
}

CureFormula::CureFormula(const InputParameters & parameters) :
  // You must call the constructor of the base class first
  Kernel(parameters),



_v_var(coupledValue("v")),// _v_var WILL HOLD TEMPERATURE VALUE

    _T_id(coupled("v")) //this is the coupled variable TEMPERATURE identifier which will be used in offDiagQPJacobian

{}

Real CureFormula::computeQpResidual()
{
  // velocity * _grad_u[_qp] is actually doing a dot product
  /* To be copied to QpJacobian */
Real a=600.0;
  Real k01=a*exp(20.2);
 // cout<<"k01 value is "<<k01<<endl;

  Real k02=a*exp(9.5);
  Real Eatt1=64735.0;
  Real Eatt2=31561.0;
  Real mm=1.0;
  Real nn=1.07;
  Real R=8.314;
  Real k1=k01*exp(-Eatt1/(R*(_v_var[_qp]+273.0)));

  Real k2=k02*exp(-Eatt2/(R*(_v_var[_qp]+273.0)));
    /* To be copied to QpJacobian */


  Real Pc=(1.0/60.0)*(k1+k2*pow(_u[_qp],mm))*pow((1.0-_u[_qp]),nn);
  //cout<<"Pc value is "<<Pc<<endl;
  return  -_test[_i][_qp]*Pc; // changed this, note - sign
//return -_test[_i][_qp]*400.0;
}

Real CureFormula::computeQpJacobian()
{
  // the partial derivative of _grad_u is just _grad_phi[_j]
    /* To be copied to QpJacobian */
Real a=600.0;
  Real k01=a*exp(20.2);
 // cout<<"k01 value is "<<k01<<endl;

  Real k02=a*exp(9.5);
  Real Eatt1=64735.0;
  Real Eatt2=31561.0;
  Real mm=1.0;
  Real nn=1.07;
  Real R=8.314;
  Real k1=k01*exp(-Eatt1/(R*(_v_var[_qp]+273.0)));

  Real k2=k02*exp(-Eatt2/(R*(_v_var[_qp]+273.0)));
    /* To be copied to QpJacobian */

  Real Pc=(1.0/60.0)*(k1+k2*pow(_u[_qp],mm))*pow((1.0-_u[_qp]),nn);
  Real dPc = (1.0/60.0)*(k2*mm*pow(_u[_qp],mm-1.0)*pow((1.0-_u[_qp]),nn)-(k1+k2*pow(_u[_qp],mm))*nn*pow((1.0-_u[_qp]),nn-1.0));
    //cout<<"dPcbydalpha value is "<<dPc<<endl;

  return -_test[_i][_qp]*dPc*_phi[_j][_qp];
  //return 0.0;
  //note that chain rule is used not in dPc calculation, but used in return statement by multiplying with _phi[_j][_qp]
  //I changed this, and note the - sign

}

//off diagonal Jacobian needs to be added !!!!!!!!!!
//
Real CureFormula::computeQpOffDiagJacobian(unsigned int jvar)
{
    cout<<"offDiagonalJacobian executed"<<endl;
    if (jvar==_T_id)
    {
/* To be copied to QpJacobian */
Real a=600.0;
  Real k01=a*exp(20.2);
 // cout<<"k01 value is "<<k01<<endl;

  Real k02=a*exp(9.5);
  Real Eatt1=64735.0;
  Real Eatt2=31561.0;
  Real mm=1.0;
  Real nn=1.07;
  Real R=8.314;
  Real k1=k01*exp(-Eatt1/(R*(_v_var[_qp]+273.0)));

  Real k2=k02*exp(-Eatt2/(R*(_v_var[_qp]+273.0)));
    /* To be copied to QpJacobian */


    Real dk1bydT=k1*(Eatt1/R)*pow((1.0/(_v_var[_qp]+273.0)),2.0);
    Real dk2bydT=k2*(Eatt2/R)*pow((1.0/(_v_var[_qp]+273.0)),2.0);

    Real dPcbydT=(1.0/60.0)*(dk1bydT+dk2bydT*pow(_u[_qp],mm))*pow((1.0-_u[_qp]),nn);
    cout<<"dPcbydT value is "<<dPcbydT<<endl;

    return -_test[_i][_qp]*dPcbydT*_phi[_j][_qp];
     //return 0.0;
     //_v_var[_qp] from residual becomes _phi[_j][_qp] in the off diagonal jacobian
    }

    return 0.0;
}
