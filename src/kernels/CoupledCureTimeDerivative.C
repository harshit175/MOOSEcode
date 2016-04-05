/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
    #include "CoupledCureTimeDerivative.h"
    using namespace std;
    template<>
    InputParameters validParams<CoupledCureTimeDerivative>()
   {
      InputParameters params = validParams<Kernel>();
      params.addClassDescription("Cure Time derivative Kernel that acts on a coupled variable");
      params.addRequiredCoupledVar("v", "Coupled variable");

      // Add a required parameter.  If this isn't provided in the input file MOOSE will error.
      // params.addRequiredParam<Real>("density", "Density of polymer");

      // Add a required parameter.  If this isn't provided in the input file MOOSE will error.
      // params.addRequiredParam<Real>("Hr", "Heat of reaction");

      return params;
    }

    CoupledCureTimeDerivative::CoupledCureTimeDerivative(const InputParameters & parameters) :
        Kernel(parameters),
_density(NULL),
	_Hr(NULL),
	_v_dot(coupledDot("v")),
	 _dv_dot(coupledDotDu("v")),
        _v_var(coupled("v"))
    {}

void
CoupledCureTimeDerivative::initialSetup()
{



    _density = &getMaterialProperty<Real>("density");
    _Hr = &getMaterialProperty<Real>("Hr");

}

    Real
    CoupledCureTimeDerivative::computeQpResidual()
    {
        //cout<<"Hr value is"<<(*_Hr)[_qp]<<endl;
      return -(*_density)[_qp]*(*_Hr)[_qp]*_test[_i][_qp] * _v_dot[_qp]; //note the - sign, i changed the sign
    }

    Real
    CoupledCureTimeDerivative::computeQpJacobian()
    {
      return 0.0;
    }

    Real
    CoupledCureTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
    {
      if (jvar == _v_var)
        return -(*_density)[_qp]*(*_Hr)[_qp]*_test[_i][_qp] * _phi[_j][_qp] * _dv_dot[_qp];

      return 0.0;
    }
