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

#include "CureFormulaS2.h"

#include <cmath>

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<CureFormulaS2>()
{
  InputParameters params = validParams<Kernel>();
params.addRequiredCoupledVar("v", "Coupled variable");
  return params;
}

CureFormulaS2::CureFormulaS2(const InputParameters & parameters) :
  // You must call the constructor of the base class first
  Kernel(parameters),


  _T_id(coupled("v")), //this is the coupled variable TEMPERATURE identifier which will be used in offDiagQPJacobian

_v_var(coupledValue("v"))// _v_var WILL HOLD TEMPERATURE VALUE
{}

Real CureFormulaS2::computeQpResidual()
{
  // velocity * _grad_u[_qp] is actually doing a dot product

  /* To be copied to QpJacobian */
  /*
  Real k01=exp(20.2);
  Real k02=exp(9.5);
  Real Eatt1=64735;
  Real Eatt2=31561;
  Real mm=1;
  Real nn=1.07;
  Real R=8.314;
  Real k1=k01*exp(-Eatt1/(R*(_v_var[_qp]+273)));
  Real k2=k02*exp(-Eatt2/(R*(_v_var[_qp]+273))); */
    /* To be copied to QpJacobian */

  /*
  Real Pc=(1/60)*(k1+k2*pow(_u[_qp],mm))*pow((1-_u[_qp]),nn);
  return  -_test[_i][_qp]*Pc; // changed this, note - sign
  */
  return  -_test[_i][_qp]*_v_var[_qp]*_u[_qp]*pow((1-_u[_qp]),3); //i changed this
    //note that _v_var[_qp] that represents linear temperature dependence has been added
}

Real CureFormulaS2::computeQpJacobian()
{
  // the partial derivative of _grad_u is just _grad_phi[_j]
    /* To be copied to QpJacobian */
  /*
  Real k01=exp(20.2);
  Real k02=exp(9.5);
  Real Eatt1=64735;
  Real Eatt2=31561;
  Real mm=1;
  Real nn=1.07;
  Real R=8.314;
  Real k1=k01*exp(-Eatt1/(R*(_v_var[_qp]+273)));
  Real k2=k02*exp(-Eatt2/(R*(_v_var[_qp]+273)));


  Real Pc=(1/60)*(k1+k2*pow(_u[_qp],mm))*pow((1-_u[_qp]),nn);
  Real dPc = (1/60)*(k2*mm*pow(_u[_qp],mm-1)*pow((1-_u[_qp]),nn)-(1/60)*(k1+k2*pow(_u[_qp],mm))*nn*pow((1-_u[_qp]),nn-1));

  return -_test[_i][_qp]*dPc*_phi[_j][_qp];
  */

  //note that chain rule is used not in dPc calculation, but used in return statement
  //I changed this, and note the - sign
  return -_test[_i][_qp]*_v_var[_qp]*(pow((1-_u[_qp]),3) +3*_u[_qp]*pow((1-_u[_qp]),2))*_phi[_j][_qp];
}

//off diagonal Jacobian needs to be added if two way coupling is there !!!!!!!!!!

Real CureFormulaS2::computeQpOffDiagJacobian(unsigned int jvar)
{
    if (jvar==_T_id)
    {
     return  -_test[_i][_qp]*_phi[_j][_qp]*_u[_qp]*pow((1-_u[_qp]),3);
     //_v_var[_qp] from residual becomes _phi[_j][_qp] in the off diagonal jacobian
    }

    return 0.0;
}
