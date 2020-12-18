//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DualRealOps.h"
#include "SymmetricIsotropicRankFourTensorImplementation.h"

template class SymmetricIsotropicRankFourTensorTempl<Real>;
template class SymmetricIsotropicRankFourTensorTempl<DualReal>;

namespace MathUtils
{
template <>
void
mooseSetToZero<SymmetricIsotropicRankFourTensor>(SymmetricIsotropicRankFourTensor & v)
{
  v.zero();
}
template <>
void
mooseSetToZero<ADSymmetricIsotropicRankFourTensor>(ADSymmetricIsotropicRankFourTensor & v)
{
  v.zero();
}
}

#define RankTwoTensorMultInstantiate(TemplateClass)                                                \
  template RankTwoTensor SymmetricIsotropicRankFourTensor::operator*(                              \
      const TemplateClass<Real> & a) const;                                                        \
  template ADRankTwoTensor ADSymmetricIsotropicRankFourTensor::operator*(                          \
      const TemplateClass<Real> & a) const;                                                        \
  template ADRankTwoTensor SymmetricIsotropicRankFourTensor::operator*(                            \
      const TemplateClass<DualReal> & a) const;                                                    \
  template ADRankTwoTensor ADSymmetricIsotropicRankFourTensor::operator*(                          \
      const TemplateClass<DualReal> & a) const

RankTwoTensorMultInstantiate(RankTwoTensorTempl);
RankTwoTensorMultInstantiate(TensorValue);
RankTwoTensorMultInstantiate(TypeTensor);

template SymmetricIsotropicRankFourTensor
SymmetricIsotropicRankFourTensor::operator+(const SymmetricIsotropicRankFourTensor & a) const;
template ADSymmetricIsotropicRankFourTensor
ADSymmetricIsotropicRankFourTensor::operator+(const SymmetricIsotropicRankFourTensor & a) const;
template ADSymmetricIsotropicRankFourTensor
SymmetricIsotropicRankFourTensor::operator+(const ADSymmetricIsotropicRankFourTensor & a) const;
template ADSymmetricIsotropicRankFourTensor
ADSymmetricIsotropicRankFourTensor::operator+(const ADSymmetricIsotropicRankFourTensor & a) const;

template SymmetricIsotropicRankFourTensor
SymmetricIsotropicRankFourTensor::operator-(const SymmetricIsotropicRankFourTensor & a) const;
template ADSymmetricIsotropicRankFourTensor
ADSymmetricIsotropicRankFourTensor::operator-(const SymmetricIsotropicRankFourTensor & a) const;
template ADSymmetricIsotropicRankFourTensor
SymmetricIsotropicRankFourTensor::operator-(const ADSymmetricIsotropicRankFourTensor & a) const;
template ADSymmetricIsotropicRankFourTensor
ADSymmetricIsotropicRankFourTensor::operator-(const ADSymmetricIsotropicRankFourTensor & a) const;

template SymmetricIsotropicRankFourTensor
    SymmetricIsotropicRankFourTensor::operator*(const SymmetricIsotropicRankFourTensor & a) const;
template ADSymmetricIsotropicRankFourTensor
    ADSymmetricIsotropicRankFourTensor::operator*(const SymmetricIsotropicRankFourTensor & a) const;
template ADSymmetricIsotropicRankFourTensor
    SymmetricIsotropicRankFourTensor::operator*(const ADSymmetricIsotropicRankFourTensor & a) const;
template ADSymmetricIsotropicRankFourTensor ADSymmetricIsotropicRankFourTensor::operator*(
    const ADSymmetricIsotropicRankFourTensor & a) const;
