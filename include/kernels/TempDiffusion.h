#ifndef TEMPDIFF_H
#define TEMPDIFF_H

#include "Diffusion.h"

class TempDiffusion;

template<>
InputParameters validParams<TempDiffusion>();


class TempDiffusion : public Diffusion
{
    public:
        TempDiffusion(const InputParameters & parameters);

    protected:
        virtual Real computeQpResidual();
        virtual Real computeQpJacobian();
	virtual void initialSetup(); //function added, to enable multi domain conductivity definition
    private:
        //Real _k; //removed
	const MaterialProperty<Real> * _k; //added,to enable multi domain conductivity definition
};

#endif  // TEMPDIFF_H
