//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolidMechanicsPropertyRegistry.h"

namespace SolidMechanics
{
Rank
operator|(Rank a, Rank b)
{
  return static_cast<Rank>(static_cast<int>(a) | static_cast<int>(b));
}

Type
operator|(Type a, Type b)
{
  return static_cast<Type>(static_cast<int>(a) | static_cast<int>(b));
}

Symmetry
operator|(Symmetry a, Symmetry b)
{
  return static_cast<Symmetry>(static_cast<int>(a) | static_cast<int>(b));
}

PropertyRegistry &
PropertyRegistry::get()
{
  static PropertyRegistry property_registry_singleton;
  return property_registry_singleton;
}

char
PropertyRegistry::add(const PropertyRegistryEntry & entry)
{
  get()._data.push_back(entry);
  return 0;
}

std::vector<PropertyRegistryEntry>
PropertyRegistry::data()
{
  return get()._data;
}

PropertyQuery
PropertyRegistry::query()
{
  return PropertyQuery();
}

PropertyQuery
PropertyQuery::filter(const PropertyRegistryEntry & e)
{
  return this->name(e.name).alias(e.alias).rank(e.rank).symmetry(e.symmetry).type(e.type);
}

PropertyQuery
PropertyQuery::name(const std::string & n)
{
  auto new_query = *this;
  if (!n.empty())
    new_query._name = n;
  return new_query;
}

PropertyQuery
PropertyQuery::alias(const std::string & a)
{
  auto new_query = *this;
  if (!a.empty())
    new_query._alias = a;
  return new_query;
}

PropertyQuery
PropertyQuery::rank(const Rank r)
{
  auto new_query = *this;
  if (r != Rank::ANY)
    new_query._rank = r;
  return new_query;
}

PropertyQuery
PropertyQuery::symmetry(const Symmetry s)
{
  auto new_query = *this;
  if (s != Symmetry::ANY)
    new_query._symmetry = s;
  return new_query;
}

PropertyQuery
PropertyQuery::type(const Type t)
{
  auto new_query = *this;
  if (t != Type::ANY)
    new_query._type = t;
  return new_query;
}

std::vector<PropertyRegistryEntry>
PropertyQuery::get() const
{
  std::vector<PropertyRegistryEntry> result;

  for (const auto & candidate : PropertyRegistry::data())
  {
    if (!_name.empty() && candidate.name != _name)
      continue;
    if (!_alias.empty() && candidate.alias != _alias)
      continue;
    if ((_rank | candidate.rank) != _rank)
      continue;
    if ((_symmetry | candidate.symmetry) != _symmetry)
      continue;
    if ((_type | candidate.type) != _type)
      continue;

    result.push_back(candidate);
  }

  return result;
}
}
