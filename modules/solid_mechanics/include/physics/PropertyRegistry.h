//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <vector>

namespace SolidMechanics
{
enum class Rank : int
{
  ZERO = 0,
  ONE = 1,
  TWO = 2,
  THREE = 4,
  FOUR = 8,
  ANY = -1
};

enum class Symmetry : int
{
  SYMMETRIC = 0,
  SKEW_SYMMETRIC = 1,
  POSSIBLY_SYMMETRIC = 2,
  NONE = 4,
  ANY = -1
};

enum class Type : int
{
  STRESS_LIKE = 0,
  STRAIN_LIKE = 1,
  ANY = -1
};

Rank operator|(Rank a, Rank b);
Symmetry operator|(Symmetry a, Symmetry b);
Type operator|(Type a, Type b);

struct PropertyRegistry
{
  /// The material property name declared inside the material object
  const std::string name;

  /// The alias used in the action for output purposes
  const std::string alias;

  /// Rank
  const Rank rank;

  /// Symmetry
  const Symmetry symmetry;

  /// Stress- or strain-like
  const Type type;
};

class QueryResult
{
public:
  QueryResult(const std::vector<PropertyRegistry> & candidates) : _candidates(candidates) {}

  QueryResult name(const std::string &);
  QueryResult alias(const std::string &);
  QueryResult rank(const Rank);
  QueryResult symmetry(const Symmetry);
  QueryResult type(const Type);

  std::vector<PropertyRegistry> get() const;

protected:
  /// @{ Query criteria
  std::string _name = "";
  std::string _alias = "";
  Rank _rank = Rank::ANY;
  Symmetry _symmetry = Symmetry::ANY;
  Type _type = Type::ANY;
  /// @}

private:
  const std::vector<PropertyRegistry> & _candidates;
};

QueryResult query(const std::vector<PropertyRegistry> & candidates);
}
