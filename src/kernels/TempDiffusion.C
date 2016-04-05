#include "TempDiffusion.h"
using namespace std;
template<>
InputParameters validParams<TempDiffusion>()
{
      InputParameters params = validParams<Diffusion>();
      // params.addRequiredParam<Real>("thermalconductivity", "Thermal Conductivity"); //removed, as now it is added via the constructor
      return params;

}

TempDiffusion::TempDiffusion(const InputParameters & parameters) :
       Diffusion(parameters), _k(NULL)
           {}

void TempDiffusion::initialSetup()
{
  _k = &getMaterialProperty<Real>("TConductivity");
}

Real TempDiffusion::computeQpResidual()
{
    //cout<<"thermal conductivity is"<<(*_k)[_qp]<<endl;
  return (*_k)[_qp]*Diffusion::computeQpResidual();
}



Real TempDiffusion::computeQpJacobian()
{
  return (*_k)[_qp]*Diffusion::computeQpJacobian();
}
