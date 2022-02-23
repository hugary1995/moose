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
  FactorizedSymmetricRankTwoTensorTempl(const RankTwoTensorTempl<T> & A) : _A(A)
  {
    A.symmetricEigenvaluesEigenvectors(_eigvals, _eigvecs);
  }

  /// Construct from RankTwoTensorTempl<T> if the factorization is known
  FactorizedSymmetricRankTwoTensorTempl(const RankTwoTensorTempl<T> & A,
                                        const std::vector<T> & eigvals,
                                        const RankTwoTensorTempl<T> & eigvecs)
    : _A(A), _eigvals(eigvals), _eigvecs(eigvecs)
  {
#ifdef DEBUG
    validate();
#endif
  }

  /// Construct from the factorization
  FactorizedSymmetricRankTwoTensorTempl(const std::vector<T> & eigvals,
                                        const RankTwoTensorTempl<T> & eigvecs)
    : _eigvals(eigvals), _eigvecs(eigvecs)
  {
    RankTwoTensorTempl<T> D(eigvals);
    _A = eigvecs * D * eigvecs.transpose();
  }

  // @{ Indexing
  const T & operator()(const unsigned int i, const unsigned int j) const { return _A(i, j); }
  T operator()(const unsigned int i, const unsigned int j) { return _A(i, j); }
  // @}

  // @{ getters
  const RankTwoTensorTempl<T> & get() const { return _A; }
  RankTwoTensorTempl<T> get() { return _A; }
  const std::vector<T> & eigvals() const { return _eigvals; }
  std::vector<T> eigvals() { return _eigvals; }
  const RankTwoTensorTempl<T> & eigvecs() const { return _eigvecs; }
  RankTwoTensorTempl<T> eigvecs() { return _eigvecs; }
  // @}

  void print(std::ostream & stm = Moose::out) const;

  /// Test if the factorization is correct, and if _A is still symmetric.
  void validate() const;

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

  /// natural log of _A
  FactorizedSymmetricRankTwoTensorTempl<T> log() const;

  /// exponentiated _A
  FactorizedSymmetricRankTwoTensorTempl<T> exp() const;

  /// _A raised to a power
  template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type = 0>
  FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype>
  pow(const T2 & p) const;

  /// sqrt of _A
  FactorizedSymmetricRankTwoTensorTempl<T> sqrt() const;

  /// cbrt of _A
  FactorizedSymmetricRankTwoTensorTempl<T> cbrt() const;

private:
  // The underlying un-factorized RankTwoTensorTempl<T>
  RankTwoTensorTempl<T> _A;

  // The eigen values of _A;
  std::vector<T> _eigvals;

  // The eigen vectors of _A;
  RankTwoTensorTempl<T> _eigvecs;
};

template <typename T>
void
FactorizedSymmetricRankTwoTensorTempl<T>::print(std::ostream & stm) const
{
  this->get().print(stm);
}

template <typename T>
void
FactorizedSymmetricRankTwoTensorTempl<T>::validate() const
{
  RankTwoTensorTempl<T> error = _A - _A.transpose();
  if (!MooseUtils::absoluteFuzzyEqual(error.norm(), 0))
    mooseError("The tensor is not symmetric.");

  RankTwoTensorTempl<T> D(_eigvals);
  RankTwoTensorTempl<T> A = _eigvecs * D * _eigvecs.transpose();
  RankTwoTensorTempl<T> error = A - _A;
  if (!MooseUtils::absoluteFuzzyEqual(error.norm(), 0))
    mooseError("Internal error: The factorization is wrong.");
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
FactorizedSymmetricRankTwoTensorTempl<T>::rotated(const RankTwoTensorTempl<T> & R) const
{
  FactorizedSymmetricRankTwoTensorTempl<T> Ar(_A.rotated(R), _eigvals, R * _eigvecs);
  return Ar;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
FactorizedSymmetricRankTwoTensorTempl<T>::transpose() const
{
  FactorizedSymmetricRankTwoTensorTempl<T> At(_A.transpose(), _eigvals, _eigvecs);
  return At;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> &
FactorizedSymmetricRankTwoTensorTempl<T>::operator=(
    const FactorizedSymmetricRankTwoTensorTempl<T> & A)
{
  _A = A._A;
  _eigvals = A._eigvals;
  _eigvecs = A._eigvecs;
  return *this;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> &
FactorizedSymmetricRankTwoTensorTempl<T>::operator=(const RankTwoTensorTempl<T> & A)
{
  _A = A;
  A.symmetricEigenvaluesEigenvectors(_eigvals, _eigvecs);
  return *this;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> &
FactorizedSymmetricRankTwoTensorTempl<T>::operator*=(const T & a)
{
  _A *= a;
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
  A._A *= a;
  for (auto & eigval : A._eigvals)
    eigval *= a;
  return A;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T> &
FactorizedSymmetricRankTwoTensorTempl<T>::operator/=(const T & a)
{
  _A /= a;
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
  A._A /= a;
  for (auto & eigval : A._eigvals)
    eigval /= a;
  return A;
}

template <typename T>
bool
FactorizedSymmetricRankTwoTensorTempl<T>::operator==(
    const FactorizedSymmetricRankTwoTensorTempl<T> & A) const
{
  for (auto i : make_range(LIBMESH_DIM))
    for (auto j : make_range(LIBMESH_DIM))
      if (!MooseUtils::absoluteFuzzyEqual((*this)(i, j), A(i, j)))
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
  _A.addIa(a);
  for (auto & eigval : _eigvals)
    eigval += a;
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
FactorizedSymmetricRankTwoTensorTempl<T>::log() const
{
  return FactorizedSymmetricRankTwoTensorTempl<T>(
      {std::log(_eigvals[0]), std::log(_eigvals[1]), std::log(_eigvals[2])}, _eigvecs);
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
FactorizedSymmetricRankTwoTensorTempl<T>::exp() const
{
  return FactorizedSymmetricRankTwoTensorTempl<T>(
      {std::exp(_eigvals[0]), std::exp(_eigvals[1]), std::exp(_eigvals[2])}, _eigvecs);
}

template <typename T>
template <typename T2, typename std::enable_if<ScalarTraits<T2>::value, int>::type>
FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype>
FactorizedSymmetricRankTwoTensorTempl<T>::pow(const T2 & p) const
{
  return FactorizedSymmetricRankTwoTensorTempl<typename CompareTypes<T, T2>::supertype>(
      {std::pow(_eigvals[0], p), std::pow(_eigvals[1], p), std::pow(_eigvals[2], p)}, _eigvecs);
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
FactorizedSymmetricRankTwoTensorTempl<T>::sqrt() const
{
  return FactorizedSymmetricRankTwoTensorTempl<T>(
      {std::sqrt(_eigvals[0]), std::sqrt(_eigvals[1]), std::sqrt(_eigvals[2])}, _eigvecs);
}

template <typename T>
FactorizedSymmetricRankTwoTensorTempl<T>
FactorizedSymmetricRankTwoTensorTempl<T>::cbrt() const
{
  return FactorizedSymmetricRankTwoTensorTempl<T>(
      {std::cbrt(_eigvals[0]), std::cbrt(_eigvals[1]), std::cbrt(_eigvals[2])}, _eigvecs);
}

typedef FactorizedSymmetricRankTwoTensorTempl<Real> FactorizedSymmetricRankTwoTensor;
typedef FactorizedSymmetricRankTwoTensorTempl<ADReal> ADFactorizedSymmetricRankTwoTensor;
