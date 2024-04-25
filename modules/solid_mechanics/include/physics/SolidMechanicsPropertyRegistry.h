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

#define registerSolidMechanicsProperty(name, rank, symmetry, type)                                 \
  static char dummyvar_for_registering_solid_mechancis_prop_##name =                               \
      SolidMechanics::PropertyRegistry::add({#name, #name, rank, symmetry, type})

#define registerSolidMechanicsPropertyAlias(name, alias, rank, symmetry, type)                     \
  static char dummyvar_for_registering_solid_mechancis_prop_##name =                               \
      SolidMechanics::PropertyRegistry::add({#name, #alias, rank, symmetry, type})

namespace SolidMechanics
{
// Forward decl
class PropertyRegistryEntry;
class PropertyRegistry;
class PropertyQuery;

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

struct PropertyRegistryEntry
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

class PropertyRegistry
{
public:
  static PropertyRegistry & get();

  static char add(const PropertyRegistryEntry & entry);

  static std::vector<PropertyRegistryEntry> data();

  static PropertyQuery query();

private:
  PropertyRegistry() {}

  std::vector<PropertyRegistryEntry> _data;
};

class PropertyQuery
{
public:
  PropertyQuery filter(const PropertyRegistryEntry &);
  PropertyQuery name(const std::string &);
  PropertyQuery alias(const std::string &);
  PropertyQuery rank(const Rank);
  PropertyQuery symmetry(const Symmetry);
  PropertyQuery type(const Type);

  std::vector<PropertyRegistryEntry> get() const;

  friend class PropertyRegistry;

protected:
  /// @{ Query criteria
  std::string _name = "";
  std::string _alias = "";
  Rank _rank = Rank::ANY;
  Symmetry _symmetry = Symmetry::ANY;
  Type _type = Type::ANY;
  /// @}

private:
  PropertyQuery() {}
};
}
