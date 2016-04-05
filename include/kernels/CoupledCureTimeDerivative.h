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

#ifndef COUPLEDCURETIMEDERIVATIVE_H
#define COUPLEDCURETIMEDERIVATIVE_H

   #include "Kernel.h"

   // Forward Declaration
   class CoupledCureTimeDerivative;

   template<>
   InputParameters validParams<CoupledCureTimeDerivative>();

   class CoupledCureTimeDerivative : public Kernel
   {
   public:
   CoupledCureTimeDerivative(const InputParameters & parameters);

   protected:
   virtual Real computeQpResidual();
   virtual Real computeQpJacobian();
   virtual Real computeQpOffDiagJacobian(unsigned int jvar);
   virtual void initialSetup(); //function added
   // Real _density; //removed
   // Real _Hr; //removed

   private:  
   const MaterialProperty<Real> * _density; //added
   const MaterialProperty<Real> * _Hr; //added
   const VariableValue & _v_dot;
   const VariableValue & _dv_dot;
   unsigned int _v_var;
   
    };

   #endif //COUPLEDCURETIMEDERIVATIVE_H
