//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeLagrangianLinearElasticStress.h"
#include "VectorizedMaterial.h"

class ComputeLagrangianLinearElasticStressVectorized : public ComputeLagrangianLinearElasticStress
{
public:
  static InputParameters validParams()
  {
    InputParameters params = ComputeLagrangianLinearElasticStress::validParams();
    params.addRequiredParam<UserObjectName>("vectorized_material", "Vectorized material");
    return params;
  }

  ComputeLagrangianLinearElasticStressVectorized(const InputParameters & parameters)
    : ComputeLagrangianLinearElasticStress(parameters),
      _vec_mat(getUserObject<VectorizedMaterial>("vectorized_material"))
  {
  }

protected:
  virtual void computeQpSmallStress() override
  {
    if (!_vec_mat.ready())
    {
      _small_stress[_qp].zero();
      _small_jacobian[_qp].zero();
      return;
    }

    _small_stress[_qp] = _vec_mat.getRankTwoTensor("small_stress", _current_elem->id(), _qp);
    _small_jacobian[_qp] = _vec_mat.getRankFourTensor("small_jacobian", _current_elem->id(), _qp);
  }

  const VectorizedMaterial & _vec_mat;
};
