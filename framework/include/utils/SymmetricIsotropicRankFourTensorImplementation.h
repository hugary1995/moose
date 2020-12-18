//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SymmetricIsotropicRankFourTensor.h"

// MOOSE includes
#include "RankTwoTensor.h"
#include "MooseEnum.h"
#include "MooseException.h"
#include "MooseUtils.h"

#include "libmesh/utility.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

// C++ includes
#include <iomanip>
#include <ostream>

namespace MathUtils
{
template <>
void mooseSetToZero<SymmetricIsotropicRankFourTensorTempl<Real>>(
    SymmetricIsotropicRankFourTensorTempl<Real> & v);
template <>
void mooseSetToZero<SymmetricIsotropicRankFourTensorTempl<DualReal>>(
    SymmetricIsotropicRankFourTensorTempl<DualReal> & v);
}

template <typename T>
MooseEnum
SymmetricIsotropicRankFourTensorTempl<T>::fillMethodEnum()
{
  return MooseEnum("antisymmetric symmetric9 symmetric21 general_isotropic symmetric_isotropic "
                   "symmetric_isotropic_E_nu antisymmetric_isotropic axisymmetric_rz general "
                   "principal");
}

template <typename T>
SymmetricIsotropicRankFourTensorTempl<T>::SymmetricIsotropicRankFourTensorTempl()
{
  mooseAssert(
      N == 3,
      "SymmetricIsotropicRankFourTensorTempl<T> is currently only tested for 3 dimensions.");

  _vals[0] = 0.;
  _vals[1] = 0.;
}

template <typename T>
SymmetricIsotropicRankFourTensorTempl<T>::SymmetricIsotropicRankFourTensorTempl(
    const InitMethod init)
{
  switch (init)
  {
    case initNone:
      break;

    case initIdentity:
      _vals[0] = 1.;
      _vals[1] = 1.;
      break;

    case initIdentityFour:
      _vals[0] = 1.;
      _vals[1] = 1.;
      break;

    case initIdentitySymmetricFour:
      _vals[0] = 1.;
      _vals[1] = 1.;
      break;

    default:
      mooseError("Unknown SymmetricIsotropicRankFourTensorTempl<T> initialization pattern.");
  }
}

template <typename T>
SymmetricIsotropicRankFourTensorTempl<T>::SymmetricIsotropicRankFourTensorTempl(
    const std::vector<T> & input, FillMethod fill_method)
{
  fillFromInputVector(input, fill_method);
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::zero()
{
  _vals[0] = 0.;
  _vals[1] = 0.;
}

template <typename T>
SymmetricIsotropicRankFourTensorTempl<T> &
SymmetricIsotropicRankFourTensorTempl<T>::operator=(
    const SymmetricIsotropicRankFourTensorTempl<T> & a)
{
  _vals[0] = a._vals[0];
  _vals[1] = a._vals[1];
  return *this;
}

template <typename T>
template <template <typename> class Tensor, typename T2>
auto SymmetricIsotropicRankFourTensorTempl<T>::operator*(const Tensor<T2> & b) const ->
    typename std::enable_if<TwoTensorMultTraits<Tensor, T2>::value,
                            RankTwoTensorTempl<decltype(T() * T2())>>::type
{
  typedef decltype(T() * T2()) ValueType;
  RankTwoTensorTempl<ValueType> result = _vals[1] * b;
  result.addIa((_vals[0] - _vals[1]) / 3. * b.tr());
  return result;
}

template <typename T>
SymmetricIsotropicRankFourTensorTempl<T> &
SymmetricIsotropicRankFourTensorTempl<T>::operator*=(const T & a)
{
  _vals[0] *= a;
  _vals[1] *= a;
  return *this;
}

template <typename T>
SymmetricIsotropicRankFourTensorTempl<T> &
SymmetricIsotropicRankFourTensorTempl<T>::operator/=(const T & a)
{
  _vals[0] /= a;
  _vals[1] /= a;
  return *this;
}

template <typename T>
SymmetricIsotropicRankFourTensorTempl<T> &
SymmetricIsotropicRankFourTensorTempl<T>::operator+=(
    const SymmetricIsotropicRankFourTensorTempl<T> & a)
{
  _vals[0] += a._vals[0];
  _vals[1] += a._vals[1];
  return *this;
}

template <typename T>
template <typename T2>
auto
SymmetricIsotropicRankFourTensorTempl<T>::operator+(
    const SymmetricIsotropicRankFourTensorTempl<T2> & b) const
    -> SymmetricIsotropicRankFourTensorTempl<decltype(T() + T2())>
{
  SymmetricIsotropicRankFourTensorTempl<decltype(T() + T2())> result;
  result._vals[0] = _vals[0] + b._vals[0];
  result._vals[1] = _vals[1] + b._vals[1];
  return result;
}

template <typename T>
SymmetricIsotropicRankFourTensorTempl<T> &
SymmetricIsotropicRankFourTensorTempl<T>::operator-=(
    const SymmetricIsotropicRankFourTensorTempl<T> & a)
{
  _vals[0] -= a._vals[0];
  _vals[1] -= a._vals[1];
  return *this;
}

template <typename T>
template <typename T2>
auto
SymmetricIsotropicRankFourTensorTempl<T>::operator-(
    const SymmetricIsotropicRankFourTensorTempl<T2> & b) const
    -> SymmetricIsotropicRankFourTensorTempl<decltype(T() - T2())>
{
  SymmetricIsotropicRankFourTensorTempl<decltype(T() - T2())> result;
  result._vals[0] = _vals[0] - b._vals[0];
  result._vals[1] = _vals[1] - b._vals[1];
  return result;
}

template <typename T>
SymmetricIsotropicRankFourTensorTempl<T>
SymmetricIsotropicRankFourTensorTempl<T>::operator-() const
{
  SymmetricIsotropicRankFourTensorTempl<T> result;
  result._vals[0] = -_vals[0];
  result._vals[1] = -_vals[1];
  return result;
}

template <typename T>
template <typename T2>
auto SymmetricIsotropicRankFourTensorTempl<T>::operator*(
    const SymmetricIsotropicRankFourTensorTempl<T2> & b) const
    -> SymmetricIsotropicRankFourTensorTempl<decltype(T() * T2())>
{
  typedef decltype(T() * T2()) ValueType;
  SymmetricIsotropicRankFourTensorTempl<ValueType> result;
  result._vals[0] = _vals[0] * b._vals[0];
  result._vals[1] = _vals[1] * b._vals[1];
  return result;
}

template <typename T>
T
SymmetricIsotropicRankFourTensorTempl<T>::L2norm() const
{
  return _vals[0] + _vals[1] * std::sqrt(5);
}

template <typename T>
SymmetricIsotropicRankFourTensorTempl<T>
SymmetricIsotropicRankFourTensorTempl<T>::invSymm() const
{
  SymmetricIsotropicRankFourTensorTempl<T> result;
  result._vals[0] = 1. / _vals[0];
  result._vals[1] = 1. / _vals[1];
  return result;
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::print(std::ostream & stm) const
{
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
    {
      stm << "i = " << i << " j = " << j << '\n';
      for (unsigned int k = 0; k < N; ++k)
      {
        for (unsigned int l = 0; l < N; ++l)
          stm << std::setw(15) << (*this)(i, j, k, l) << " ";

        stm << '\n';
      }
    }
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::surfaceFillFromInputVector(
    const std::vector<T> & /*input*/)
{
  mooseError("not yet implemented");
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillFromInputVector(const std::vector<T> & input,
                                                              FillMethod fill_method)
{

  switch (fill_method)
  {
    case antisymmetric:
      fillAntisymmetricFromInputVector(input);
      break;
    case symmetric9:
      fillSymmetric9FromInputVector(input);
      break;
    case symmetric21:
      fillSymmetric21FromInputVector(input);
      break;
    case general_isotropic:
      fillGeneralIsotropicFromInputVector(input);
      break;
    case symmetric_isotropic:
      fillSymmetricIsotropicFromInputVector(input);
      break;
    case symmetric_isotropic_E_nu:
      fillSymmetricIsotropicEandNuFromInputVector(input);
      break;
    case antisymmetric_isotropic:
      fillAntisymmetricIsotropicFromInputVector(input);
      break;
    case axisymmetric_rz:
      fillAxisymmetricRZFromInputVector(input);
      break;
    case general:
      fillGeneralFromInputVector(input);
      break;
    case principal:
      fillPrincipalFromInputVector(input);
      break;
    default:
      mooseError("fillFromInputVector called with unknown fill_method of ", fill_method);
  }
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillSymmetric9FromInputVector(
    const std::vector<T> & /*input*/)
{
  mooseError("not yet implemented");
}
template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillSymmetric21FromInputVector(
    const std::vector<T> & /*input*/)
{
  mooseError("not yet implemented");
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillAntisymmetricFromInputVector(
    const std::vector<T> & /*input*/)
{
  mooseError("not yet implemented");
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillGeneralIsotropicFromInputVector(
    const std::vector<T> & /*input*/)
{
  mooseError("not yet implemented");
}

template <typename T>
void SymmetricIsotropicRankFourTensorTempl<T>::fillGeneralIsotropic(T /*i0*/, T /*i1*/, T /*i2*/)
{
  mooseError("not yet implemented");
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillAntisymmetricIsotropicFromInputVector(
    const std::vector<T> & /*input*/)
{
  mooseError("not yet implemented");
}

template <typename T>
void SymmetricIsotropicRankFourTensorTempl<T>::fillAntisymmetricIsotropic(T /*i0*/)
{
  mooseError("not yet implemented");
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillSymmetricIsotropicFromInputVector(
    const std::vector<T> & input)
{
  fillSymmetricIsotropic(input[0], input[1]);
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillSymmetricIsotropic(T K, T G)
{
  _vals[0] = 3. * K;
  _vals[1] = 2. * G;
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillSymmetricIsotropicEandNuFromInputVector(
    const std::vector<T> & input)
{
  fillSymmetricIsotropicEandNu(input[0], input[1]);
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillSymmetricIsotropicEandNu(T E, T nu)
{
  // Calculate lambda and the shear modulus from the given young's modulus and poisson's ratio
  const T & K = E / 3. / (1. - 2. * nu);
  const T & G = E / (2.0 * (1.0 + nu));

  fillSymmetricIsotropic(K, G);
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillAxisymmetricRZFromInputVector(
    const std::vector<T> & /*input*/)
{
  mooseError("not yet implemented");
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillGeneralFromInputVector(
    const std::vector<T> & /*input*/)
{
  mooseError("not yet implemented");
}

template <typename T>
void
SymmetricIsotropicRankFourTensorTempl<T>::fillPrincipalFromInputVector(
    const std::vector<T> & /*input*/)
{
  mooseError("not yet implemented");
}

template <typename T>
RankTwoTensorTempl<T>
SymmetricIsotropicRankFourTensorTempl<T>::innerProductTranspose(
    const RankTwoTensorTempl<T> & /*b*/) const
{
  mooseError("not yet implemented");
  RankTwoTensorTempl<T> result;
  return result;
}

template <typename T>
T
SymmetricIsotropicRankFourTensorTempl<T>::sum3x3() const
{
  mooseError("not yet implemented");
  return 0.;
}

template <typename T>
VectorValue<T>
SymmetricIsotropicRankFourTensorTempl<T>::sum3x1() const
{
  mooseError("not yet implemented");
  VectorValue<T> a(3);
  return a;
}

template <typename T>
bool
SymmetricIsotropicRankFourTensorTempl<T>::isSymmetric() const
{
  return true;
}

template <typename T>
bool
SymmetricIsotropicRankFourTensorTempl<T>::isIsotropic() const
{
  return true;
}
