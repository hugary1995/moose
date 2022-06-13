//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralUserObject.h"
#include "DependencyResolverInterface.h"
#include "Conversion.h"

// TIMPI includes
#include "timpi/communicator.h"
#include "timpi/parallel_sync.h"

/* This class gathers material properties requested by a material object to form a vector */
class VectorizedMaterialBase : public ElementUserObject, public DependencyResolverInterface
{
public:
  static InputParameters validParams();

  VectorizedMaterialBase(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject &) override;
  virtual const std::set<std::string> & getRequestedItems() override { return _reqs; }
  virtual const std::set<std::string> & getSuppliedItems() override { return _sups; }

  // The public accessor -- only searches for the material property in the _output vector.
  template <typename T>
  T get(const std::string name, unsigned int elem, unsigned int qp) const
  {
    return get<T>(_output, name, elem, qp);
  }

  virtual bool ready() const { return _ready; }

protected:
  // This is the entry point for GPU calls.
  // You could safely assume that _input and _input_old vectors are all gathered on rank 0.
  // Use fill and get methods below to load and store material properties at a specific element/qp.
  virtual void GPUCalls() = 0;

  unsigned int dof(unsigned int elem, unsigned int qp) const
  {
    return elem * _material_data->nQPoints() + qp;
  }

  void allocateVector(std::map<std::string, std::vector<Real>> & container,
                      const std::set<std::string> & prop_names);

  void fillVector(std::map<std::string, std::vector<Real>> & container,
                  const std::set<std::string> & prop_names,
                  const MaterialProperties & props);

  template <typename T>
  inline void fill(std::map<std::string, std::vector<Real>> & container,
                   const std::string name,
                   unsigned int elem,
                   unsigned int qp,
                   const T & data);

  template <typename T>
  inline T get(const std::map<std::string, std::vector<Real>> & container,
               const std::string name,
               unsigned int elem,
               unsigned int qp) const;

  MaterialBase * _mat;
  std::map<MaterialPropertyName, std::string> _type;

  std::map<std::string, unsigned int> _prop_name_id_map;
  std::map<unsigned int, std::string> _prop_id_name_map;

  mutable std::set<std::string> _reqs;
  std::set<std::string> _sups;

  std::set<std::string> _input_prop_names;
  std::set<std::string> _input_old_prop_names;
  std::set<std::string> _output_prop_names;

  std::map<std::string, std::vector<Real>> _input;
  std::map<std::string, std::vector<Real>> _input_old;
  std::map<std::string, std::vector<Real>> _output;

private:
  void
  serialize(std::vector<Real> & vec,
            const std::unordered_map<processor_id_type, std::vector<dof_id_type>> & query_dofs);

  void
  deserialize(std::vector<Real> & vec,
              const std::unordered_map<processor_id_type, std::vector<dof_id_type>> & query_dofs);

  template <typename T>
  inline void allocateVectorHelper(std::map<std::string, std::vector<Real>> & container,
                                   const std::string prop_name);

  template <typename T>
  const MooseArray<T> & getProp(const MaterialProperties & props, const std::string prop_name) const
  {
    return static_cast<MaterialPropertyBase<T, false> *>(props[_prop_name_id_map.at(prop_name)])
        ->get();
  }

  bool _ready;
};

template <>
inline void
VectorizedMaterialBase::allocateVectorHelper<Real>(
    std::map<std::string, std::vector<Real>> & container, const std::string prop_name)
{
  container[prop_name].resize(_mesh.maxElemId() * _material_data->nQPoints());
}

template <>
inline void
VectorizedMaterialBase::allocateVectorHelper<RankTwoTensor>(
    std::map<std::string, std::vector<Real>> & container, const std::string prop_name)
{
  for (auto i : make_range(3))
    for (auto j : make_range(3))
      container[prop_name + "_" + Moose::stringify(i) + Moose::stringify(j)].resize(
          _mesh.maxElemId() * _material_data->nQPoints());
}

template <>
inline void
VectorizedMaterialBase::allocateVectorHelper<RankFourTensor>(
    std::map<std::string, std::vector<Real>> & container, const std::string prop_name)
{
  for (auto i : make_range(3))
    for (auto j : make_range(3))
      for (auto k : make_range(3))
        for (auto l : make_range(3))
          container[prop_name + "_" + Moose::stringify(i) + Moose::stringify(j) +
                    Moose::stringify(k) + Moose::stringify(l)]
              .resize(_mesh.maxElemId() * _material_data->nQPoints());
}

template <>
inline void
VectorizedMaterialBase::fill<Real>(std::map<std::string, std::vector<Real>> & container,
                                   const std::string name,
                                   unsigned int elem,
                                   unsigned int qp,
                                   const Real & data)
{
  container[name][dof(elem, qp)] = data;
}

template <>
inline void
VectorizedMaterialBase::fill<RankTwoTensor>(std::map<std::string, std::vector<Real>> & container,
                                            const std::string name,
                                            unsigned int elem,
                                            unsigned int qp,
                                            const RankTwoTensor & data)
{
  for (auto i : make_range(3))
    for (auto j : make_range(3))
      container[name + "_" + Moose::stringify(i) + Moose::stringify(j)][dof(elem, qp)] = data(i, j);
}

template <>
inline void
VectorizedMaterialBase::fill<RankFourTensor>(std::map<std::string, std::vector<Real>> & container,
                                             const std::string name,
                                             unsigned int elem,
                                             unsigned int qp,
                                             const RankFourTensor & data)
{
  for (auto i : make_range(3))
    for (auto j : make_range(3))
      for (auto k : make_range(3))
        for (auto l : make_range(3))
          container[name + "_" + Moose::stringify(i) + Moose::stringify(j) + Moose::stringify(k) +
                    Moose::stringify(l)][dof(elem, qp)] = data(i, j, k, l);
}

template <>
inline Real
VectorizedMaterialBase::get<Real>(const std::map<std::string, std::vector<Real>> & container,
                                  const std::string name,
                                  unsigned int elem,
                                  unsigned int qp) const
{
  return container.at(name)[dof(elem, qp)];
}

template <>
inline RankTwoTensor
VectorizedMaterialBase::get<RankTwoTensor>(
    const std::map<std::string, std::vector<Real>> & container,
    const std::string name,
    unsigned int elem,
    unsigned int qp) const
{
  RankTwoTensor out;
  for (auto i : make_range(3))
    for (auto j : make_range(3))
      out(i, j) =
          container.at(name + "_" + Moose::stringify(i) + Moose::stringify(j))[dof(elem, qp)];
  return out;
}

template <>
inline RankFourTensor
VectorizedMaterialBase::get<RankFourTensor>(
    const std::map<std::string, std::vector<Real>> & container,
    const std::string name,
    unsigned int elem,
    unsigned int qp) const
{
  RankFourTensor out;
  for (auto i : make_range(3))
    for (auto j : make_range(3))
      for (auto k : make_range(3))
        for (auto l : make_range(3))
          out(i, j, k, l) = container.at(name + "_" + Moose::stringify(i) + Moose::stringify(j) +
                                         Moose::stringify(k) + Moose::stringify(l))[dof(elem, qp)];
  return out;
}
