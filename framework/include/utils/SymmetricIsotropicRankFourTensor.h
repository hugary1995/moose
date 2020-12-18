//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Moose.h"
#include "ADRankTwoTensorForward.h"
#include "ADRankFourTensorForward.h"

#include "libmesh/libmesh.h"

#include "metaphysicl/raw_type.h"

using libMesh::Real;
namespace libMesh
{
template <typename>
class TensorValue;
template <typename>
class TypeTensor;
template <typename>
class VectorValue;
}

// Forward declarations
class MooseEnum;

namespace MathUtils
{
template <typename T>
void mooseSetToZero(T & v);

/**
 * Helper function template specialization to set an object to zero.
 * Needed by DerivativeMaterialInterface
 */
template <>
void mooseSetToZero<SymmetricIsotropicRankFourTensor>(SymmetricIsotropicRankFourTensor & v);
template <>
void mooseSetToZero<ADSymmetricIsotropicRankFourTensor>(ADSymmetricIsotropicRankFourTensor & v);
}

/**
 * SymmetricIsotropicRankFourTensorTempl is the efficient and economical implementation of a
 * symmetric isotropic RankFourTensor
 */
template <typename T>
class SymmetricIsotropicRankFourTensorTempl
{
public:
  /// Initialization method
  enum InitMethod
  {
    initNone,
    initIdentity,
    initIdentityFour,
    initIdentitySymmetricFour
  };

  /**
   * To fill up the 81 entries in the 4th-order tensor, fillFromInputVector
   * is called with one of the following fill_methods.
   * See the fill*FromInputVector functions for more details
   */
  enum FillMethod
  {
    antisymmetric,
    symmetric9,
    symmetric21,
    general_isotropic,
    symmetric_isotropic,
    symmetric_isotropic_E_nu,
    antisymmetric_isotropic,
    axisymmetric_rz,
    general,
    principal
  };

  template <template <typename> class Tensor, typename Scalar>
  struct TwoTensorMultTraits
  {
    static const bool value = false;
  };
  template <typename Scalar>
  struct TwoTensorMultTraits<RankTwoTensorTempl, Scalar>
  {
    static const bool value = ScalarTraits<Scalar>::value;
  };
  template <typename Scalar>
  struct TwoTensorMultTraits<TensorValue, Scalar>
  {
    static const bool value = ScalarTraits<Scalar>::value;
  };
  template <typename Scalar>
  struct TwoTensorMultTraits<TypeTensor, Scalar>
  {
    static const bool value = ScalarTraits<Scalar>::value;
  };

  /// Default constructor; fills to zero
  SymmetricIsotropicRankFourTensorTempl();

  /// Select specific initialization pattern
  SymmetricIsotropicRankFourTensorTempl(const InitMethod);

  /// Fill from vector
  SymmetricIsotropicRankFourTensorTempl(const std::vector<T> &, FillMethod);

  /// Copy assignment operator must be defined if used
  SymmetricIsotropicRankFourTensorTempl(const SymmetricIsotropicRankFourTensorTempl<T> & a) =
      default;

  /**
   * Copy constructor
   */
  template <typename T2>
  SymmetricIsotropicRankFourTensorTempl(const SymmetricIsotropicRankFourTensorTempl<T2> & copy);

  // Named constructors
  static SymmetricIsotropicRankFourTensorTempl<T> Identity()
  {
    return SymmetricIsotropicRankFourTensorTempl<T>(initIdentity);
  }
  static SymmetricIsotropicRankFourTensorTempl<T> IdentityFour()
  {
    return SymmetricIsotropicRankFourTensorTempl<T>(initIdentityFour);
  };

  /**
   * Gets the value for the index specified.  Takes index = 0,1,2
   * used for const
   */
  inline T operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const
  {
    const T I = (i == j) * (k == l) / 3.;
    const T J = ((i == k) * (j == l) + (i == l) * (j == k)) / 2. - I;
    return _vals[0] * I + _vals[1] * J;
  }

  /**
   * Gets the principal components
   */
  inline T & operator()(unsigned int i) { return _vals[i]; }
  inline T operator()(unsigned int i) const { return _vals[i]; }

  /// Zeros out the tensor.
  void zero();

  /// Print the rank four tensor
  void print(std::ostream & stm = Moose::out) const;

  /// copies values from a into this tensor
  SymmetricIsotropicRankFourTensorTempl<T> &
  operator=(const SymmetricIsotropicRankFourTensorTempl<T> & a);

  /**
   * Assignment-from-scalar operator.  Used only to zero out the tensor.
   *
   * \returns A reference to *this.
   */
  template <typename Scalar>
  typename boostcopy::enable_if_c<ScalarTraits<Scalar>::value,
                                  SymmetricIsotropicRankFourTensorTempl &>::type
  operator=(const Scalar & libmesh_dbg_var(p))
  {
    libmesh_assert_equal_to(p, Scalar(0));
    this->zero();
    return *this;
  }

  /// C_ijkl*a_kl
  template <template <typename> class Tensor, typename T2>
  auto operator*(const Tensor<T2> & a) const ->
      typename std::enable_if<TwoTensorMultTraits<Tensor, T2>::value,
                              RankTwoTensorTempl<decltype(T() * T2())>>::type;

  /// C_ijkl*a
  template <typename T2>
  auto operator*(const T2 & a) const ->
      typename std::enable_if<ScalarTraits<T2>::value,
                              SymmetricIsotropicRankFourTensorTempl<decltype(T() * T2())>>::type;

  /// C_ijkl *= a
  SymmetricIsotropicRankFourTensorTempl<T> & operator*=(const T & a);

  /// C_ijkl/a
  template <typename T2>
  auto operator/(const T2 & a) const ->
      typename std::enable_if<ScalarTraits<T2>::value,
                              SymmetricIsotropicRankFourTensorTempl<decltype(T() / T2())>>::type;

  /// C_ijkl /= a  for all i, j, k, l
  SymmetricIsotropicRankFourTensorTempl<T> & operator/=(const T & a);

  /// C_ijkl += a_ijkl  for all i, j, k, l
  SymmetricIsotropicRankFourTensorTempl<T> &
  operator+=(const SymmetricIsotropicRankFourTensorTempl<T> & a);

  /// C_ijkl + a_ijkl
  template <typename T2>
  auto operator+(const SymmetricIsotropicRankFourTensorTempl<T2> & a) const
      -> SymmetricIsotropicRankFourTensorTempl<decltype(T() + T2())>;

  /// C_ijkl -= a_ijkl
  SymmetricIsotropicRankFourTensorTempl<T> &
  operator-=(const SymmetricIsotropicRankFourTensorTempl<T> & a);

  /// C_ijkl - a_ijkl
  template <typename T2>
  auto operator-(const SymmetricIsotropicRankFourTensorTempl<T2> & a) const
      -> SymmetricIsotropicRankFourTensorTempl<decltype(T() - T2())>;

  /// -C_ijkl
  SymmetricIsotropicRankFourTensorTempl<T> operator-() const;

  /// C_ijpq*a_pqkl
  template <typename T2>
  auto operator*(const SymmetricIsotropicRankFourTensorTempl<T2> & a) const
      -> SymmetricIsotropicRankFourTensorTempl<decltype(T() * T2())>;

  /// sqrt(C_ijkl*C_ijkl)
  T L2norm() const;

  /// inv(C)
  SymmetricIsotropicRankFourTensorTempl<T> invSymm() const;

  /// rotate an isotropic tensor
  void rotate(const TypeTensor<T> &) {}

  /// transpose a symmetric tensor
  SymmetricIsotropicRankFourTensorTempl<T> transposeMajor() const { return *this; }

  /**
   * Fills the tensor entries ignoring the last dimension (ie, C_ijkl=0 if any of i, j, k, or l =
   * 3).
   * Fill method depends on size of input
   * Input size = 2.  Then C_1111 = C_2222 = input[0], and C_1122 = input[1], and C_1212 = (input[0]
   * - input[1])/2,
   *                  and C_ijkl = C_jikl = C_ijlk = C_klij, and C_1211 = C_1222 = 0.
   * Input size = 9.  Then C_1111 = input[0], C_1112 = input[1], C_1122 = input[3],
   *                       C_1212 = input[4], C_1222 = input[5], C_1211 = input[6]
   *                       C_2211 = input[7], C_2212 = input[8], C_2222 = input[9]
   *                       and C_ijkl = C_jikl = C_ijlk
   */
  void surfaceFillFromInputVector(const std::vector<T> & input);

  /// Static method for use in validParams for getting the "fill_method"
  static MooseEnum fillMethodEnum();

  /**
   * fillFromInputVector takes some number of inputs to fill
   * the Rank-4 tensor.
   * @param input the numbers that will be placed in the tensor
   * @param fill_method this can be:
   *             antisymmetric (use fillAntisymmetricFromInputVector)
   *             symmetric9 (use fillSymmetric9FromInputVector)
   *             symmetric21 (use fillSymmetric21FromInputVector)
   *             general_isotropic (use fillGeneralIsotropicFrominputVector)
   *             symmetric_isotropic (use fillSymmetricIsotropicFromInputVector)
   *             antisymmetric_isotropic (use fillAntisymmetricIsotropicFromInputVector)
   *             axisymmetric_rz (use fillAxisymmetricRZFromInputVector)
   *             general (use fillGeneralFromInputVector)
   *             principal (use fillPrincipalFromInputVector)
   */
  void fillFromInputVector(const std::vector<T> & input, FillMethod fill_method);

  ///@{ Vector-less fill API functions. See docs of the corresponding ...FromInputVector methods
  void fillGeneralIsotropic(T i0, T i1, T i2);
  void fillAntisymmetricIsotropic(T i0);
  void fillSymmetricIsotropic(T i0, T i1);
  void fillSymmetricIsotropicEandNu(T E, T nu);
  ///@}

  /// Inner product of the major transposed tensor with a rank two tensor
  RankTwoTensorTempl<T> innerProductTranspose(const RankTwoTensorTempl<T> &) const;

  /// Calculates the sum of Ciijj for i and j varying from 0 to 2
  T sum3x3() const;

  /// Calculates the vector a[i] = sum over j Ciijj for i and j varying from 0 to 2
  VectorValue<T> sum3x1() const;

  /// checks if the tensor is symmetric
  bool isSymmetric() const;

  /// checks if the tensor is isotropic
  bool isIsotropic() const;

protected:
  /// Dimensionality of rank-four tensor
  static constexpr unsigned int N = LIBMESH_DIM;
  static constexpr unsigned int N2 = N * N;
  static constexpr unsigned int N3 = N * N * N;
  static constexpr unsigned int N4 = N * N * N * N;

  /**
   * C = _vals[0] * I + _vals[1] * J
   * where I is the volumetric fourth order identity projector, J is the deviatoric projector
   * I = 1/3 * delta_ij * delta_kl
   * J = 1/2 * (delta_ik * delta_jl + delta_il * delta_jk) - I
   * I and J are both projectors, i.e. I * I = I, J * J = J
   * I and J are orthogonal, i.e. I * J = 0
   */
  T _vals[2];

  /**
   * fillSymmetric9FromInputVector takes 9 inputs to fill in
   * the Rank-4 tensor with the appropriate crystal symmetries maintained. I.e., C_ijkl = C_klij,
   * C_ijkl = C_ijlk, C_ijkl = C_jikl
   * @param input is:
   *                C1111 C1122 C1133 C2222 C2233 C3333 C2323 C1313 C1212
   *                In the isotropic case this is (la is first Lame constant, mu is second (shear)
   * Lame constant)
   *                la+2mu la la la+2mu la la+2mu mu mu mu
   */
  void fillSymmetric9FromInputVector(const std::vector<T> & input);

  /**
   * fillSymmetric21FromInputVector takes either 21 inputs to fill in
   * the Rank-4 tensor with the appropriate crystal symmetries maintained. I.e., C_ijkl = C_klij,
   * C_ijkl = C_ijlk, C_ijkl = C_jikl
   * @param input is
   *                C1111 C1122 C1133 C1123 C1113 C1112 C2222 C2233 C2223 C2213 C2212 C3333 C3323
   * C3313 C3312 C2323 C2313 C2312 C1313 C1312 C1212
   */
  void fillSymmetric21FromInputVector(const std::vector<T> & input);

  /**
   * fillAntisymmetricFromInputVector takes 6 inputs to fill the
   * the antisymmetric Rank-4 tensor with the appropriate symmetries maintained.
   * I.e., B_ijkl = -B_jikl = -B_ijlk = B_klij
   * @param input this is B1212, B1213, B1223, B1313, B1323, B2323.
   */
  void fillAntisymmetricFromInputVector(const std::vector<T> & input);

  /**
   * fillGeneralIsotropicFromInputVector takes 3 inputs to fill the
   * Rank-4 tensor with symmetries C_ijkl = C_klij, and isotropy, ie
   * C_ijkl = la*de_ij*de_kl + mu*(de_ik*de_jl + de_il*de_jk) + a*ep_ijm*ep_klm
   * where la is the first Lame modulus, mu is the second (shear) Lame modulus,
   * and a is the antisymmetric shear modulus, and ep is the permutation tensor
   * @param input this is la, mu, a in the above formula
   */
  void fillGeneralIsotropicFromInputVector(const std::vector<T> & input);

  /**
   * fillAntisymmetricIsotropicFromInputVector takes 1 input to fill the
   * the antisymmetric Rank-4 tensor with the appropriate symmetries maintained.
   * I.e., C_ijkl = a * ep_ijm * ep_klm, where epsilon is the permutation tensor (and sum on m)
   * @param input this is a in the above formula
   */
  void fillAntisymmetricIsotropicFromInputVector(const std::vector<T> & input);

  /**
   * fillSymmetricIsotropicFromInputVector takes 2 inputs to fill the
   * the symmetric Rank-4 tensor with the appropriate symmetries maintained.
   * C_ijkl = lambda*de_ij*de_kl + mu*(de_ik*de_jl + de_il*de_jk)
   * where lambda is the first Lame modulus, mu is the second (shear) Lame modulus,
   * @param input this is lambda and mu in the above formula
   */
  void fillSymmetricIsotropicFromInputVector(const std::vector<T> & input);

  /**
   * fillSymmetricIsotropicEandNuFromInputVector is a variation of the
   * fillSymmetricIsotropicFromInputVector which takes as inputs the
   * more commonly used Young's modulus (E) and Poisson's ratio (nu)
   * constants to fill the isotropic elasticity tensor. Using well-known formulas,
   * E and nu are used to calculate lambda and mu and then the vector is passed
   * to fillSymmetricIsotropicFromInputVector.
   * @param input Young's modulus (E) and Poisson's ratio (nu)
   */
  void fillSymmetricIsotropicEandNuFromInputVector(const std::vector<T> & input);

  /**
   * fillAxisymmetricRZFromInputVector takes 5 inputs to fill the axisymmetric
   * Rank-4 tensor with the appropriate symmetries maintatined for use with
   * axisymmetric problems using coord_type = RZ.
   * I.e. C1111 = C2222, C1133 = C2233, C2323 = C3131 and C1212 = 0.5*(C1111-C1122)
   * @param input this is C1111, C1122, C1133, C3333, C2323.
   */
  void fillAxisymmetricRZFromInputVector(const std::vector<T> & input);

  /**
   * fillGeneralFromInputVector takes 81 inputs to fill the Rank-4 tensor
   * No symmetries are explicitly maintained
   * @param input  C(i,j,k,l) = input[i*N*N*N + j*N*N + k*N + l]
   */
  void fillGeneralFromInputVector(const std::vector<T> & input);

  /**
   * fillPrincipalFromInputVector takes 9 inputs to fill a Rank-4 tensor
   * C1111 = input0
   * C1122 = input1
   * C1133 = input2
   * C2211 = input3
   * C2222 = input4
   * C2233 = input5
   * C3311 = input6
   * C3322 = input7
   * C3333 = input8
   * with all other components being zero
   */

  void fillPrincipalFromInputVector(const std::vector<T> & input);
  template <class T2>
  friend void dataStore(std::ostream &, SymmetricIsotropicRankFourTensorTempl<T2> &, void *);

  template <class T2>
  friend void dataLoad(std::istream &, SymmetricIsotropicRankFourTensorTempl<T2> &, void *);

  template <typename T2>
  friend class RankTwoTensorTempl;
  template <typename T2>
  friend class SymmetricIsotropicRankFourTensorTempl;
  template <typename T2>
  friend class RankThreeTensorTempl;
};

namespace MetaPhysicL
{
template <typename T>
struct RawType<SymmetricIsotropicRankFourTensorTempl<T>>
{
  typedef SymmetricIsotropicRankFourTensorTempl<typename RawType<T>::value_type> value_type;

  static value_type value(const SymmetricIsotropicRankFourTensorTempl<T> & in)
  {
    value_type ret;
    ret(0) = raw_value(in(0));
    ret(1) = raw_value(in(1));
    return ret;
  }
};
}

template <typename T1, typename T2>
inline auto operator*(const T1 & a, const SymmetricIsotropicRankFourTensorTempl<T2> & b) ->
    typename std::enable_if<ScalarTraits<T1>::value,
                            SymmetricIsotropicRankFourTensorTempl<decltype(T1() * T2())>>::type
{
  return b * a;
}

template <typename T>
template <typename T2>
SymmetricIsotropicRankFourTensorTempl<T>::SymmetricIsotropicRankFourTensorTempl(
    const SymmetricIsotropicRankFourTensorTempl<T2> & copy)
{
  _vals[0] = copy._vals[0];
  _vals[1] = copy._vals[1];
}

template <typename T>
template <typename T2>
auto SymmetricIsotropicRankFourTensorTempl<T>::operator*(const T2 & b) const ->
    typename std::enable_if<ScalarTraits<T2>::value,
                            SymmetricIsotropicRankFourTensorTempl<decltype(T() * T2())>>::type
{
  typedef decltype(T() * T2()) ValueType;
  SymmetricIsotropicRankFourTensorTempl<ValueType> result;

  result._vals[0] = _vals[0] * b;
  result._vals[1] = _vals[1] * b;

  return result;
}

template <typename T>
template <typename T2>
auto
SymmetricIsotropicRankFourTensorTempl<T>::operator/(const T2 & b) const ->
    typename std::enable_if<ScalarTraits<T2>::value,
                            SymmetricIsotropicRankFourTensorTempl<decltype(T() / T2())>>::type
{
  SymmetricIsotropicRankFourTensorTempl<decltype(T() / T2())> result;
  result._vals[0] = _vals[0] / b;
  result._vals[1] = _vals[1] / b;
  return result;
}
