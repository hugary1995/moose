//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolidMechanicsMaterialProperty.h"

using namespace SolidMechanicsMaterialProperty;

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

QueryResult
QueryResult::name(const std::string & n)
{
  auto new_query = *this;
  new_query._name = n;
  return new_query;
}

QueryResult
QueryResult::alias(const std::string & a)
{
  auto new_query = *this;
  new_query._alias = a;
  return new_query;
}

QueryResult
QueryResult::rank(const Rank r)
{
  auto new_query = *this;
  new_query._rank = r;
  return new_query;
}

QueryResult
QueryResult::symmetry(const Symmetry s)
{
  auto new_query = *this;
  new_query._symmetry = s;
  return new_query;
}

QueryResult
QueryResult::type(const Type t)
{
  auto new_query = *this;
  new_query._type = t;
  return new_query;
}

std::vector<Registry>
QueryResult::get() const
{
  std::vector<Registry> result;

  for (const auto & candidate : _candidates)
  {
    if (!_name.empty() && candidate.name != _name)
      continue;
    if (!_alias.empty() && candidate.alias != _alias)
      continue;
    if (_rank | candidate.rank != _rank)
      continue;
    if (_symmetry | candidate.symmetry != _symmetry)
      continue;
    if (_type | candidate.type != _type)
      continue;
    result.push_back(candidate);
  }

  return result;
}

QueryResult
query(const std::vector<Registry> & candidates)
{
  return QueryResult(candidates);
}
