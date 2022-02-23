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
class FactorizedSymmetricRankTwoTensorTempl;

namespace MathUtils
{
/// natural log of a factorized RankTwoTensor
template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> log(const FactorizedSymmetricRankTwoTensorTempl<T> &);

/// exponentiated a factorized RankTwoTensor
template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> exp(const FactorizedSymmetricRankTwoTensorTempl<T> &);

/// a factorized RankTwoTensor raised to a power
template <typename T, typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type = 0>
FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype>
pow(const FactorizedSymmetricRankTwoTensorTempl<T> &, const T2 & p);

/// sqrt of a factorized RankTwoTensor
template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> sqrt(const FactorizedSymmetricRankTwoTensorTempl<T> &);

/// cbrt of a factorized RankTwoTensor
template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> cbrt(const FactorizedSymmetricRankTwoTensorTempl<T> &);
} // end namespace MathUtils

/**
 * FactorizedSymmetricRankTwoTensorTempl is designed to perform the spectral decomposition of an
 * underlying RankTwoTensorTempl and reuse its bases for future operations if possible.
 *
 * Only operations that reuses the known factorization are provided. Otherwise, you will need to
 * first retrieve the underlying RankTwoTensorTempl to perform the operation.
 */
template <typename T>
class FactorizedSymmetricRankTwoTensorTempl
{
public:
  /// No default constructor
  FactorizedSymmetricRankTwoTensorTempl() = delete;

  /// Copy constructor
  FactorizedSymmetricRankTwoTensorTempl(const FactorizedSymmetricRankTwoTensorTempl<T> & A) =
      default;

  /// Construct from RankTwoTensorTempl<T> if the factorization isn't known a priori
  FactorizedSymmetricRankTwoTensorTempl(const RankTwoTensorTempl<T> & A)
  {
#ifdef DEBUG
    RankTwoTensorTempl<T> error = A - A.transpose();
    if (!MooseUtils::absoluteFuzzyEqual(error.norm(), 0))
      mooseError("The tensor is not symmetric.");
#endif
    A.symmetricEigenvaluesEigenvectors(_eigvals, _eigvecs);
  }

  /// Construct from RankTwoTensorTempl<T> if the factorization is known
  FactorizedSymmetricRankTwoTensorTempl(const std::vector<T> & eigvals,
                                        const RankTwoTensorTempl<T> & eigvecs)
    : _eigvals(eigvals), _eigvecs(eigvecs)
  {
  }

  // @{ getters
  RankTwoTensorTempl<T> get() const { return assemble(); }
  RankTwoTensorTempl<T> get() { return assemble(); }
  const std::vector<T> & eigvals() const { return _eigvals; }
  std::vector<T> eigvals() { return _eigvals; }
  const RankTwoTensorTempl<T> & eigvecs() const { return _eigvecs; }
  RankTwoTensorTempl<T> eigvecs() { return _eigvecs; }
  // @}

  void print(std::ostream & stm = Moose::out) const;

  /// Returns _A rotated by R.
  FactorizedSymmetricRankTwoTensorTempl<T> rotated(const RankTwoTensorTempl<T> & R) const;

  /// Returns the transpose of _A.
  FactorizedSymmetricRankTwoTensorTempl<T> transpose() const;

  // @{ Assignment operators
  FactorizedSymmetricRankTwoTensorTempl<T> &
  operator=(const FactorizedSymmetricRankTwoTensorTempl<T> & A);
  FactorizedSymmetricRankTwoTensorTempl<T> & operator=(const RankTwoTensorTempl<T> & A);
  // @}

  /// performs _A *= a, also updates eigen values
  FactorizedSymmetricRankTwoTensorTempl<T> & operator*=(const T & a);

  /// returns _A * a, also updates eigen values
  template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type = 0>
  FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype>
  operator*(const T2 & a) const;

  /// performs _A /= a, also updates eigen values
  FactorizedSymmetricRankTwoTensorTempl<T> & operator/=(const T & a);

  /// returns _A / a, also updates eigen values
  template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type = 0>
  FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype>
  operator/(const T2 & a) const;

  /// Defines logical equality with another RankTwoTensorTempl<T>
  bool operator==(const FactorizedSymmetricRankTwoTensorTempl<T> & A) const;

  /// inverse of _A
  FactorizedSymmetricRankTwoTensorTempl<T> inverse() const;

  /// add identity times a to _A
  void addIa(const T & a);

private:
  // Assemble the tensor from the factorization
  RankTwoTensorTempl<T> assemble() const
  {
    RankTwoTensorTempl<T> D(_eigvals);
    return _eigvecs * D * _eigvecs.transpose();
  }

  // The eigen values of _A;
  std::vector<T> _eigvals;

  // The eigen vectors of _A;
  RankTwoTensorTempl<T> _eigvecs;
};

namespace MathUtils
{
template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
log(const FactorizedSymmetricRankTwoTensorTempl<T> & A)
{
  const auto & eigvals = A.eigvals();
  return FactorizedSymmetricRankTwoTensorTempl<T>(
      {std::log(eigvals[0]), std::log(eigvals[1]), std::log(eigvals[2])}, A.eigvecs());
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
exp(const FactorizedSymmetricRankTwoTensorTempl<T> & A)
{
  const auto & eigvals = A.eigvals();
  return FactorizedSymmetricRankTwoTensorTempl<T>(
      {std::exp(eigvals[0]), std::exp(eigvals[1]), std::exp(eigvals[2])}, A.eigvecs());
}

template <typename T, typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type>
FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype>
pow(const FactorizedSymmetricRankTwoTensorTempl<T> & A, const T2 & p)
{
  const auto & eigvals = A.eigvals();
  return FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype>(
      {std::pow(eigvals[0], p), std::pow(eigvals[1], p), std::pow(eigvals[2], p)}, A.eigvecs());
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
sqrt(const FactorizedSymmetricRankTwoTensorTempl<T> & A)
{
  const auto & eigvals = A.eigvals();
  return FactorizedSymmetricRankTwoTensorTempl<T>(
      {std::sqrt(eigvals[0]), std::sqrt(eigvals[1]), std::sqrt(eigvals[2])}, A.eigvecs());
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
cbrt(const FactorizedSymmetricRankTwoTensorTempl<T> & A)
{
  const auto & eigvals = A.eigvals();
  return FactorizedSymmetricRankTwoTensorTempl<T>(
      {std::cbrt(eigvals[0]), std::cbrt(eigvals[1]), std::cbrt(eigvals[2])}, A.eigvecs());
}
} // end namespace MathUtils

template <typename T>
void
FactorizedSymmetricRankTwoTensorTempl<T>::print(std::ostream & stm) const
{
  this->get().print(stm);
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
FactorizedSymmetricRankTwoTensorTempl<T>::rotated(const RankTwoTensorTempl<T> & R) const
{
  return FactorizedSymmetricRankTwoTensorTempl<T>(_eigvals, R * _eigvecs);
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
FactorizedSymmetricRankTwoTensorTempl<T>::transpose() const
{
  return *this;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> &
FactorizedSymmetricRankTwoTensorTempl<T>::operator=(
    const FactorizedSymmetricRankTwoTensorTempl<T> & A)
{
  _eigvals = A._eigvals;
  _eigvecs = A._eigvecs;
  return *this;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> &
FactorizedSymmetricRankTwoTensorTempl<T>::operator=(const RankTwoTensorTempl<T> & A)
{
  A.symmetricEigenvaluesEigenvectors(_eigvals, _eigvecs);
  return *this;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> &
FactorizedSymmetricRankTwoTensorTempl<T>::operator*=(const T & a)
{
  for (auto & eigval : _eigvals)
    eigval *= a;
  return *this;
}

template <typename T>
template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type>
FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype>
FactorizedSymmetricRankTwoTensorTempl<T>::operator*(const T2 & a) const
{
  FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype> A = *this;
  for (auto & eigval : A._eigvals)
    eigval *= a;
  return A;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> &
FactorizedSymmetricRankTwoTensorTempl<T>::operator/=(const T & a)
{
  for (auto & eigval : _eigvals)
    eigval /= a;
  return *this;
}

template <typename T>
template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type>
FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype>
FactorizedSymmetricRankTwoTensorTempl<T>::operator/(const T2 & a) const
{
  FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype> A = *this;
  for (auto & eigval : A._eigvals)
    eigval /= a;
  return A;
}

template <typename T>
bool
FactorizedSymmetricRankTwoTensorTempl<T>::operator==(
    const FactorizedSymmetricRankTwoTensorTempl<T> & A) const
{
  RankTwoTensorTempl<T> me = get();
  RankTwoTensorTempl<T> you = A.get();
  for (auto i : make_range(LIBMESH_DIM))
    for (auto j : make_range(LIBMESH_DIM))
      if (!MooseUtils::absoluteFuzzyEqual(me(i, j), you(i, j)))
        return false;

  return true;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
FactorizedSymmetricRankTwoTensorTempl<T>::inverse() const
{
  return FactorizedSymmetricRankTwoTensorTempl<T>(
      {1 / _eigvals[0], 1 / _eigvals[1], 1 / _eigvals[2]}, _eigvecs);
}

template <typename T>
void
FactorizedSymmetricRankTwoTensorTempl<T>::addIa(const T & a)
{
  for (auto & eigval : _eigvals)
    eigval += a;
}

typedef FactorizedSymmetricRankTwoTensorTempl<Real> FactorizedSymmetricRankTwoTensor;
typedef FactorizedSymmetricRankTwoTensorTempl<ADReal> ADFactorizedSymmetricRankTwoTensor;
