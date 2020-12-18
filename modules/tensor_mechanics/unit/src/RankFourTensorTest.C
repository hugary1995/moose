//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "SymmetricIsotropicRankFourTensor.h"
#include "MooseTypes.h"
#include "ADReal.h"

#include "metaphysicl/raw_type.h"

RankFourTensor iSymmetric = RankFourTensor(RankFourTensor::initIdentitySymmetricFour);

TEST(RankFourTensor, invSymm1)
{
  // inverse check for a standard symmetric-isotropic tensor
  std::vector<Real> input(2);
  input[0] = 1;
  input[1] = 3;
  RankFourTensor a(input, RankFourTensor::symmetric_isotropic);

  EXPECT_NEAR(0, (iSymmetric - a.invSymm() * a).L2norm(), 1E-5);
}

TEST(RankFourTensor, invSymm2)
{
  // following (basically random) "a" tensor has symmetry
  // a_ijkl = a_jikl = a_ijlk
  // BUT it doesn't have a_ijkl = a_klij
  RankFourTensor a;
  a(0, 0, 0, 0) = 1;
  a(0, 0, 0, 1) = a(0, 0, 1, 0) = 2;
  a(0, 0, 0, 2) = a(0, 0, 2, 0) = 1.1;
  a(0, 0, 1, 1) = 0.5;
  a(0, 0, 1, 2) = a(0, 0, 2, 1) = -0.2;
  a(0, 0, 2, 2) = 0.8;
  a(1, 1, 0, 0) = 0.6;
  a(1, 1, 0, 1) = a(1, 1, 1, 0) = 1.3;
  a(1, 1, 0, 2) = a(1, 1, 2, 0) = 0.9;
  a(1, 1, 1, 1) = 0.6;
  a(1, 1, 1, 2) = a(1, 1, 2, 1) = -0.3;
  a(1, 1, 2, 2) = -1.1;
  a(2, 2, 0, 0) = -0.6;
  a(2, 2, 0, 1) = a(2, 2, 1, 0) = -0.1;
  a(2, 2, 0, 2) = a(2, 2, 2, 0) = -0.3;
  a(2, 2, 1, 1) = 0.5;
  a(2, 2, 1, 2) = a(2, 2, 2, 1) = 0.7;
  a(2, 2, 2, 2) = -0.9;
  a(0, 1, 0, 0) = a(1, 0, 0, 0) = 1.2;
  a(0, 1, 0, 1) = a(0, 1, 1, 0) = a(1, 0, 0, 1) = a(1, 0, 1, 0) = 0.1;
  a(0, 1, 0, 2) = a(0, 1, 2, 0) = a(1, 0, 0, 2) = a(1, 0, 2, 0) = 0.15;
  a(0, 1, 1, 1) = a(1, 0, 1, 1) = 0.24;
  a(0, 1, 1, 2) = a(0, 1, 2, 1) = a(1, 0, 1, 2) = a(1, 0, 2, 1) = -0.4;
  a(0, 1, 2, 2) = a(1, 0, 2, 2) = -0.3;
  a(0, 2, 0, 0) = a(2, 0, 0, 0) = -0.3;
  a(0, 2, 0, 1) = a(0, 2, 1, 0) = a(2, 0, 0, 1) = a(2, 0, 1, 0) = -0.4;
  a(0, 2, 0, 2) = a(0, 2, 2, 0) = a(2, 0, 0, 2) = a(2, 0, 2, 0) = 0.22;
  a(0, 2, 1, 1) = a(2, 0, 1, 1) = -0.9;
  a(0, 2, 1, 2) = a(0, 2, 2, 1) = a(2, 0, 1, 2) = a(2, 0, 2, 1) = -0.05;
  a(0, 2, 2, 2) = a(2, 0, 2, 2) = 0.32;
  a(1, 2, 0, 0) = a(2, 1, 0, 0) = -0.35;
  a(1, 2, 0, 1) = a(1, 2, 1, 0) = a(2, 1, 0, 1) = a(2, 1, 1, 0) = -1.4;
  a(1, 2, 0, 2) = a(1, 2, 2, 0) = a(2, 1, 0, 2) = a(2, 1, 2, 0) = 0.2;
  a(1, 2, 1, 1) = a(2, 1, 1, 1) = -0.91;
  a(1, 2, 1, 2) = a(1, 2, 2, 1) = a(2, 1, 1, 2) = a(2, 1, 2, 1) = 1.4;
  a(1, 2, 2, 2) = a(2, 1, 2, 2) = 0.1;

  EXPECT_NEAR(0, (iSymmetric - a.invSymm() * a).L2norm(), 1E-5);
}

TEST(RankFourTensor, ADConversion)
{
  RankFourTensor reg;
  ADRankFourTensor ad;

  ad = reg;
  reg = MetaPhysicL::raw_value(ad);

  GenericRankFourTensor<false> generic_reg;
  GenericRankFourTensor<true> generic_ad;

  generic_ad = generic_reg;
  generic_reg = MetaPhysicL::raw_value(generic_ad);
}

TEST(RankFourTensor, anisotropic21)
{

  std::vector<Real> input(21);
  for (size_t i = 0; i < input.size(); ++i)
    input[i] = i + 1.1;

  RankFourTensor a(input, RankFourTensor::symmetric21);
  RankFourTensor b;
  b(0, 0, 0, 0) = 1.1;
  b(0, 0, 0, 1) = 6.1;
  b(0, 0, 0, 2) = 5.1;
  b(0, 0, 1, 0) = 6.1;
  b(0, 0, 1, 1) = 2.1;
  b(0, 0, 1, 2) = 4.1;
  b(0, 0, 2, 0) = 5.1;
  b(0, 0, 2, 1) = 4.1;
  b(0, 0, 2, 2) = 3.1;
  b(0, 1, 0, 0) = 6.1;
  b(0, 1, 0, 1) = 21.1;
  b(0, 1, 0, 2) = 20.1;
  b(0, 1, 1, 0) = 21.1;
  b(0, 1, 1, 1) = 11.1;
  b(0, 1, 1, 2) = 18.1;
  b(0, 1, 2, 0) = 20.1;
  b(0, 1, 2, 1) = 18.1;
  b(0, 1, 2, 2) = 15.1;
  b(0, 2, 0, 0) = 5.1;
  b(0, 2, 0, 1) = 20.1;
  b(0, 2, 0, 2) = 19.1;
  b(0, 2, 1, 0) = 20.1;
  b(0, 2, 1, 1) = 10.1;
  b(0, 2, 1, 2) = 17.1;
  b(0, 2, 2, 0) = 19.1;
  b(0, 2, 2, 1) = 17.1;
  b(0, 2, 2, 2) = 14.1;
  b(1, 0, 0, 0) = 6.1;
  b(1, 0, 0, 1) = 21.1;
  b(1, 0, 0, 2) = 20.1;
  b(1, 0, 1, 0) = 21.1;
  b(1, 0, 1, 1) = 11.1;
  b(1, 0, 1, 2) = 18.1;
  b(1, 0, 2, 0) = 20.1;
  b(1, 0, 2, 1) = 18.1;
  b(1, 0, 2, 2) = 15.1;
  b(1, 1, 0, 0) = 2.1;
  b(1, 1, 0, 1) = 11.1;
  b(1, 1, 0, 2) = 10.1;
  b(1, 1, 1, 0) = 11.1;
  b(1, 1, 1, 1) = 7.1;
  b(1, 1, 1, 2) = 9.1;
  b(1, 1, 2, 0) = 10.1;
  b(1, 1, 2, 1) = 9.1;
  b(1, 1, 2, 2) = 8.1;
  b(1, 2, 0, 0) = 4.1;
  b(1, 2, 0, 1) = 18.1;
  b(1, 2, 0, 2) = 17.1;
  b(1, 2, 1, 0) = 18.1;
  b(1, 2, 1, 1) = 9.1;
  b(1, 2, 1, 2) = 16.1;
  b(1, 2, 2, 0) = 17.1;
  b(1, 2, 2, 1) = 16.1;
  b(1, 2, 2, 2) = 13.1;
  b(2, 0, 0, 0) = 5.1;
  b(2, 0, 0, 1) = 20.1;
  b(2, 0, 0, 2) = 19.1;
  b(2, 0, 1, 0) = 20.1;
  b(2, 0, 1, 1) = 10.1;
  b(2, 0, 1, 2) = 17.1;
  b(2, 0, 2, 0) = 19.1;
  b(2, 0, 2, 1) = 17.1;
  b(2, 0, 2, 2) = 14.1;
  b(2, 1, 0, 0) = 4.1;
  b(2, 1, 0, 1) = 18.1;
  b(2, 1, 0, 2) = 17.1;
  b(2, 1, 1, 0) = 18.1;
  b(2, 1, 1, 1) = 9.1;
  b(2, 1, 1, 2) = 16.1;
  b(2, 1, 2, 0) = 17.1;
  b(2, 1, 2, 1) = 16.1;
  b(2, 1, 2, 2) = 13.1;
  b(2, 2, 0, 0) = 3.1;
  b(2, 2, 0, 1) = 15.1;
  b(2, 2, 0, 2) = 14.1;
  b(2, 2, 1, 0) = 15.1;
  b(2, 2, 1, 1) = 8.1;
  b(2, 2, 1, 2) = 13.1;
  b(2, 2, 2, 0) = 14.1;
  b(2, 2, 2, 1) = 13.1;
  b(2, 2, 2, 2) = 12.1;
  EXPECT_NEAR(0, (a - b).L2norm(), 1E-5);
}

TEST(RankFourTensor, anisotropic9)
{

  std::vector<Real> input(9);
  for (size_t i = 0; i < input.size(); ++i)
    input[i] = i + 1.1;

  RankFourTensor a(input, RankFourTensor::symmetric9);
  RankFourTensor b;
  b(0, 0, 0, 0) = 1.1;
  b(0, 0, 0, 1) = 0;
  b(0, 0, 0, 2) = 0;
  b(0, 0, 1, 0) = 0;
  b(0, 0, 1, 1) = 2.1;
  b(0, 0, 1, 2) = 0;
  b(0, 0, 2, 0) = 0;
  b(0, 0, 2, 1) = 0;
  b(0, 0, 2, 2) = 3.1;
  b(0, 1, 0, 0) = 0;
  b(0, 1, 0, 1) = 9.1;
  b(0, 1, 0, 2) = 0;
  b(0, 1, 1, 0) = 9.1;
  b(0, 1, 1, 1) = 0;
  b(0, 1, 1, 2) = 0;
  b(0, 1, 2, 0) = 0;
  b(0, 1, 2, 1) = 0;
  b(0, 1, 2, 2) = 0;
  b(0, 2, 0, 0) = 0;
  b(0, 2, 0, 1) = 0;
  b(0, 2, 0, 2) = 8.1;
  b(0, 2, 1, 0) = 0;
  b(0, 2, 1, 1) = 0;
  b(0, 2, 1, 2) = 0;
  b(0, 2, 2, 0) = 8.1;
  b(0, 2, 2, 1) = 0;
  b(0, 2, 2, 2) = 0;
  b(1, 0, 0, 0) = 0;
  b(1, 0, 0, 1) = 9.1;
  b(1, 0, 0, 2) = 0;
  b(1, 0, 1, 0) = 9.1;
  b(1, 0, 1, 1) = 0;
  b(1, 0, 1, 2) = 0;
  b(1, 0, 2, 0) = 0;
  b(1, 0, 2, 1) = 0;
  b(1, 0, 2, 2) = 0;
  b(1, 1, 0, 0) = 2.1;
  b(1, 1, 0, 1) = 0;
  b(1, 1, 0, 2) = 0;
  b(1, 1, 1, 0) = 0;
  b(1, 1, 1, 1) = 4.1;
  b(1, 1, 1, 2) = 0;
  b(1, 1, 2, 0) = 0;
  b(1, 1, 2, 1) = 0;
  b(1, 1, 2, 2) = 5.1;
  b(1, 2, 0, 0) = 0;
  b(1, 2, 0, 1) = 0;
  b(1, 2, 0, 2) = 0;
  b(1, 2, 1, 0) = 0;
  b(1, 2, 1, 1) = 0;
  b(1, 2, 1, 2) = 7.1;
  b(1, 2, 2, 0) = 0;
  b(1, 2, 2, 1) = 7.1;
  b(1, 2, 2, 2) = 0;
  b(2, 0, 0, 0) = 0;
  b(2, 0, 0, 1) = 0;
  b(2, 0, 0, 2) = 8.1;
  b(2, 0, 1, 0) = 0;
  b(2, 0, 1, 1) = 0;
  b(2, 0, 1, 2) = 0;
  b(2, 0, 2, 0) = 8.1;
  b(2, 0, 2, 1) = 0;
  b(2, 0, 2, 2) = 0;
  b(2, 1, 0, 0) = 0;
  b(2, 1, 0, 1) = 0;
  b(2, 1, 0, 2) = 0;
  b(2, 1, 1, 0) = 0;
  b(2, 1, 1, 1) = 0;
  b(2, 1, 1, 2) = 7.1;
  b(2, 1, 2, 0) = 0;
  b(2, 1, 2, 1) = 7.1;
  b(2, 1, 2, 2) = 0;
  b(2, 2, 0, 0) = 3.1;
  b(2, 2, 0, 1) = 0;
  b(2, 2, 0, 2) = 0;
  b(2, 2, 1, 0) = 0;
  b(2, 2, 1, 1) = 5.1;
  b(2, 2, 1, 2) = 0;
  b(2, 2, 2, 0) = 0;
  b(2, 2, 2, 1) = 0;
  b(2, 2, 2, 2) = 6.1;
  EXPECT_NEAR(0, (a - b).L2norm(), 1E-5);
}

TEST(RankFourTensor, rotation)
{
  std::vector<Real> input(81);
  for (size_t i = 0; i < input.size(); ++i)
    input[i] = i + 1.1;

  RankFourTensor a(input, RankFourTensor::general);

  RotationTensor xrot(RotationTensor::XAXIS, 90);
  RotationTensor yrot(RotationTensor::YAXIS, 90);

  auto xyrot = xrot * yrot;
  auto yxrot = yrot * xrot;

  auto axy1 = a;
  axy1.rotate(xrot);
  axy1.rotate(yrot);

  auto axy2 = a;
  axy2.rotate(yxrot);

  EXPECT_NEAR(0, (axy1 - axy2).L2norm(), 1E-8);
}

TEST(SymmetricIsotropicRankFourTensor, indexing)
{
  Real lambda = 1.532;
  Real G = 3.132;
  Real K = lambda + 2. / 3. * G;
  std::vector<Real> input(2);

  // general rank four tensor
  input[0] = lambda;
  input[1] = G;
  RankFourTensor a(input, RankFourTensor::symmetric_isotropic);

  // symmetric isotropic rank four tensor
  input[0] = K;
  input[1] = G;
  SymmetricIsotropicRankFourTensor b(input, SymmetricIsotropicRankFourTensor::symmetric_isotropic);

  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
          EXPECT_NEAR(a(i, j, k, l), b(i, j, k, l), 1E-6);
}

TEST(SymmetricIsotropicRankFourTensor, multiply)
{
  Real lambda = 1.532;
  Real G = 3.132;
  Real K = lambda + 2. / 3. * G;
  std::vector<Real> input(2);

  // general rank four tensor
  input[0] = lambda;
  input[1] = G;
  RankFourTensor a(input, RankFourTensor::symmetric_isotropic);

  // symmetric isotropic rank four tensor
  input[0] = K;
  input[1] = G;
  SymmetricIsotropicRankFourTensor b(input, SymmetricIsotropicRankFourTensor::symmetric_isotropic);

  // a non symmetric rank two tensor
  RankTwoTensor c(1.235, 3.231, 6.452, 7.567, 3.231, 4.323, 6.432, 8.654, 9.241);
  RankTwoTensor ac = a * c;
  RankTwoTensor bc = b * c;

  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      EXPECT_NEAR(ac(i, j), bc(i, j), 1E-6);
}

TEST(SymmetricIsotropicRankFourTensor, inverse)
{
  Real lambda = 1.532;
  Real G = 3.132;
  Real K = lambda + 2. / 3. * G;
  std::vector<Real> input(2);

  // general rank four tensor
  input[0] = lambda;
  input[1] = G;
  RankFourTensor a(input, RankFourTensor::symmetric_isotropic);

  // symmetric isotropic rank four tensor
  input[0] = K;
  input[1] = G;
  SymmetricIsotropicRankFourTensor b(input, SymmetricIsotropicRankFourTensor::symmetric_isotropic);

  // a non symmetric rank two tensor
  RankFourTensor ainv = a.invSymm();
  SymmetricIsotropicRankFourTensor binv = b.invSymm();

  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
          EXPECT_NEAR(ainv(i, j, k, l), binv(i, j, k, l), 1E-6);
}
