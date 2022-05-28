//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RankTwoTensor.h"

// forward declarations
template <typename>
class FactorizedRankTwoTensorTempl;

// namespace MathUtils
// {
// /// natural log of a factorized RankTwoTensor
// template <typename T>
// FactorizedRankTwoTensorTempl<T> log(const FactorizedRankTwoTensorTempl<T> &);
//
// /// exponentiated a factorized RankTwoTensor
// template <typename T>
// FactorizedRankTwoTensorTempl<T> exp(const FactorizedRankTwoTensorTempl<T> &);
//
// /// sqrt of a factorized RankTwoTensor
// template <typename T>
// FactorizedRankTwoTensorTempl<T> sqrt(const FactorizedRankTwoTensorTempl<T> &);
//
// /// cbrt of a factorized RankTwoTensor
// template <typename T>
// FactorizedRankTwoTensorTempl<T> cbrt(const FactorizedRankTwoTensorTempl<T> &);
//
// /// a factorized RankTwoTensor raised to a power
// template <typename T, typename T2>
// FactorizedRankTwoTensorTempl<T> pow(const FactorizedRankTwoTensorTempl<T> &, const T2 &);
// } // end namespace MathUtils

/**
 * FactorizedRankTwoTensorTempl is designed to perform the spectral decomposition of an
 * underlying symmetric second order tensor and reuse its bases for future operations if possible.
 *
 * FactorizedRankTwoTensorTempl templates on the underlying second order tensor.
 * IMPORTANT: the underlying second order tensor must be symmetric. A check is only performed in
 * debug mode.
 *
 * Only operations that reuses the known factorization are provided. Otherwise, you will need to
 * first retrieve the underlying RankTwoTensorTempl to perform the operation.
 *
 * TODO? Although I am calling it factorization, it only really refers to eigenvalue decompositon at
 * this point. I am not sure if in the future we need to add other types of similarity
 * transformation.
 */
template <typename T>
class FactorizedRankTwoTensorTempl
{
public:
  /// For generic programming
  typedef typename T::value_type value_type;

  /// No default constructor
  FactorizedRankTwoTensorTempl() = delete;

  /// Copy constructor
  FactorizedRankTwoTensorTempl(const FactorizedRankTwoTensorTempl<T> & A) = default;

  /// Constructor if the factorization isn't known a priori
  FactorizedRankTwoTensorTempl(const T & A)
  {
#ifdef DEBUG
    if (!A.isSymmetric())
      mooseError("The tensor is not symmetric.");
#endif
    A.symmetricEigenvaluesEigenvectors(_eigvals, _eigvecs);
  }

  // Construct from the factorization
  // Assume that regardless of the type of T, the eigenvectors are always stored as as a full second
  // order tensor, e.g. the eigenvector matrix of a symmetric second order tensor isn't symmetric in
  // general. Hence we don't take advantage of orthonormal and/or unitary second order tensors.
  FactorizedRankTwoTensorTempl(const std::vector<typename T::value_type> & eigvals,
                               const RankTwoTensorTempl<typename T::value_type> & eigvecs)
    : _eigvals(eigvals), _eigvecs(eigvecs)
  {
  }

  // @{ getters
  template <typename T2 = RankTwoTensorTempl<typename T::value_type>>
  T2 get() const
  {
    return static_cast<T2>(assemble());
  }
  template <typename T2 = RankTwoTensorTempl<typename T::value_type>>
  T2 get()
  {
    return static_cast<T2>(assemble());
  }
  const std::vector<typename T::value_type> & eigvals() const { return _eigvals; }
  std::vector<typename T::value_type> eigvals() { return _eigvals; }
  const RankTwoTensorTempl<typename T::value_type> & eigvecs() const { return _eigvecs; }
  RankTwoTensorTempl<typename T::value_type> eigvecs() { return _eigvecs; }
  // @}

  void print(std::ostream & stm = Moose::out) const;

  // Returns _A rotated by R.
  FactorizedRankTwoTensorTempl<T>
  rotated(const RankTwoTensorTempl<typename T::value_type> & R) const;

  /// Returns the transpose of _A.
  FactorizedRankTwoTensorTempl<T> transpose() const;

  // @{ Assignment operators
  FactorizedRankTwoTensorTempl<T> & operator=(const FactorizedRankTwoTensorTempl<T> & A);
  FactorizedRankTwoTensorTempl<T> & operator=(const T & A);
  // @}

  /// performs _A *= a in place, also updates eigen values
  FactorizedRankTwoTensorTempl<T> & operator*=(const typename T::value_type & a);

  /// returns _A * a, also updates eigen values
  template <typename T2>
  FactorizedRankTwoTensorTempl<T> operator*(const T2 & a) const;

  /// performs _A /= a in place, also updates eigen values
  FactorizedRankTwoTensorTempl<T> & operator/=(const typename T::value_type & a);

  /// returns _A / a, also updates eigen values
  template <typename T2>
  FactorizedRankTwoTensorTempl<T> operator/(const T2 & a) const;

  /// Defines logical equality with another second order tensor
  bool operator==(const T & A) const;

  /// Defines logical equality with another FactorizedRankTwoTensorTempl<T>
  bool operator==(const FactorizedRankTwoTensorTempl<T> & A) const;

  /// inverse of _A
  FactorizedRankTwoTensorTempl<T> inverse() const;

  /// add identity times a to _A
  void addIa(const typename T::value_type & a);

private:
  // Assemble the tensor from the factorization
  RankTwoTensorTempl<typename T::value_type> assemble() const
  {
    RankTwoTensorTempl<typename T::value_type> D(_eigvals);
    return _eigvecs * D * _eigvecs.transpose();
  }

  // The eigen values of _A;
  std::vector<typename T::value_type> _eigvals;

  // The eigen vectors of _A;
  RankTwoTensorTempl<typename T::value_type> _eigvecs;
};

namespace MathUtils
{
#define FactorizedRankTwoTensorOperatorMapStdUnary(operator)                                       \
  template <typename T>                                                                            \
  FactorizedRankTwoTensorTempl<T> operator(const FactorizedRankTwoTensorTempl<T> & A)              \
  {                                                                                                \
    const auto & eigvals = A.eigvals();                                                            \
    return FactorizedRankTwoTensorTempl<T>(                                                        \
        {std::operator(eigvals[0]), std::operator(eigvals[1]), std::operator(eigvals[2])},         \
        A.eigvecs());                                                                              \
  }

#define FactorizedRankTwoTensorOperatorMapStdBinary(operator)                                      \
  template <typename T, typename T2>                                                               \
  FactorizedRankTwoTensorTempl<T> operator(const FactorizedRankTwoTensorTempl<T> & A,              \
                                           const T2 & b)                                           \
  {                                                                                                \
    if constexpr (ScalarTraits<T2>::value)                                                         \
    {                                                                                              \
      const auto & eigvals = A.eigvals();                                                          \
      return FactorizedRankTwoTensorTempl<T>({std::operator(eigvals[0], b),                        \
                                              std::operator(eigvals[1], b),                        \
                                              std::operator(eigvals[2], b)},                       \
                                             A.eigvecs());                                         \
    }                                                                                              \
  }

// TODO: While the macro is here, in the future we could instantiate other operator maps like
// trignometry functions.
FactorizedRankTwoTensorOperatorMapStdUnary(log);
FactorizedRankTwoTensorOperatorMapStdUnary(exp);
FactorizedRankTwoTensorOperatorMapStdUnary(sqrt);
FactorizedRankTwoTensorOperatorMapStdUnary(cbrt);
FactorizedRankTwoTensorOperatorMapStdBinary(pow);
} // end namespace MathUtils

template <typename T>
template <typename T2>
FactorizedRankTwoTensorTempl<T>
FactorizedRankTwoTensorTempl<T>::operator*(const T2 & a) const
{
  if constexpr (ScalarTraits<T2>::value)
  {
    FactorizedRankTwoTensorTempl<T> A = *this;
    for (auto & eigval : A._eigvals)
      eigval *= a;
    return A;
  }
}

template <typename T>
template <typename T2>
FactorizedRankTwoTensorTempl<T>
FactorizedRankTwoTensorTempl<T>::operator/(const T2 & a) const
{
  if constexpr (ScalarTraits<T2>::value)
  {
    FactorizedRankTwoTensorTempl<T> A = *this;
    for (auto & eigval : A._eigvals)
      eigval /= a;
    return A;
  }
}

typedef FactorizedRankTwoTensorTempl<RankTwoTensor> FactorizedRankTwoTensor;
typedef FactorizedRankTwoTensorTempl<ADRankTwoTensor> ADFactorizedRankTwoTensor;
// In the future, we could also instantiate on
// typedef FactorizedRankTwoTensorTempl<SymmtericRankTwoTensor> FactorizedSymmetricRankTwoTensor;
// typedef FactorizedRankTwoTensorTempl<ADSymmetricRankTwoTensor>
// ADFactorizedSymmetricRankTwoTensor;
