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

#include "ExampleConvection.h"
#include <cmath>

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<ExampleConvection>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<Real>("source", "Heat Source Term");
  params.addRequiredParam<Real>("A","Arrhenius constant");
  return params;
}

ExampleConvection::ExampleConvection(const InputParameters & parameters) :
  // You must call the constructor of the base class first
  Kernel(parameters),
  _source(getParam<Real>("source")),_A(getParam<Real>("A"))
{}

Real ExampleConvection::computeQpResidual()
{
  // velocity * _grad_u[_qp] is actually doing a dot product
  return _test[_i][_qp]*pow((2.71828),(-_A/_u[_qp]))*_source; //i changed this
  //return _test[_i][_qp]*_source;
}
Real ExampleConvection::computeQpJacobian()
{
  // the partial derivative of _grad_u is just _grad_phi[_j]
  // the partial derivative of _source is 0.
  return _test[_i][_qp]*(_A/pow((_u[_qp]),2))*pow((2.71828),(-_A/_u[_qp]))*_source*_phi[_j][_qp]; //i changed this
  //return 0;
}
