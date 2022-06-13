//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VectorizedMaterialBase.h"
#include "MooseMesh.h"

InputParameters
VectorizedMaterialBase::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addRequiredParam<MaterialName>("material", "The material to vectorize.");
  return params;
}

VectorizedMaterialBase::VectorizedMaterialBase(const InputParameters & parameters)
  : ElementUserObject(parameters), _ready(false)
{
  _sups.insert(name());
}

void
VectorizedMaterialBase::initialSetup()
{
  _mat = &getMaterialByName(getParam<MaterialName>("material"), /*no_warn=*/true);

  // Useful property name, id, type maps
  const auto & storage = _material_data->getMaterialPropertyStorage();
  const auto & prop_name_id_map = storage.propIDs();
  for (const auto & [name, id] : prop_name_id_map)
  {
    _prop_name_id_map.emplace(name, id);
    _prop_id_name_map.emplace(id, name);
    _type.emplace(name, _material_data->props()[id]->type());
  }

  // Figure out names of material properties to vectorize.
  const auto & reqs = _mat->getRequestedItems();
  const auto & sups = _mat->getSuppliedItems();
  const auto & dep_ids = _mat->getMatPropDependencies();
  _input_prop_names.insert(reqs.begin(), reqs.end());
  _output_prop_names.insert(sups.begin(), sups.end());
  for (const auto & dep_id : dep_ids)
    if (storage.isStatefulProp(_prop_id_name_map[dep_id]))
      _input_old_prop_names.insert(_prop_id_name_map[dep_id]);

  // The requested material properties must be reinited before the execution of this userobject.
  _reqs.insert(reqs.begin(), reqs.end());
  _material_property_dependencies.insert(dep_ids.begin(), dep_ids.end());

  // Debugging output
  _console << COLOR_CYAN << "    Input material property names: " << COLOR_DEFAULT;
  for (const auto & name : _input_prop_names)
    _console << name << " ";
  _console << std::endl;
  _console << COLOR_CYAN << "Input old material property names: " << COLOR_DEFAULT;
  for (const auto & name : _input_old_prop_names)
    _console << name << " ";
  _console << std::endl;
  _console << COLOR_CYAN << "   Output material property names: " << COLOR_DEFAULT;
  for (const auto & name : _output_prop_names)
    _console << name << " ";
  _console << std::endl;
}

void
VectorizedMaterialBase::initialize()
{
  _ready = false;

  // Components are stored separately, e.g.
  // Material property of type Real is stored as 1 vector of Real
  // Material property of type RankTwoTensor is stored as 9 vectors of Real
  // Material property of type RankFourTensor is stored as 81 vectors of Real
  allocateVector(_input, _input_prop_names);
  allocateVector(_input_old, _input_old_prop_names);
  allocateVector(_output, _output_prop_names);
}

void
VectorizedMaterialBase::execute()
{
  fillVector(_input, _input_prop_names, _material_data->props());
  fillVector(_input_old, _input_old_prop_names, _material_data->propsOld());
}

void
VectorizedMaterialBase::finalize()
{
  // We do all of the following with strong assumptions:
  // 1. There is only 1 node
  // 2. There is only 1 GPU on this node
  // Therefore, we gather all input vectors onto rank 0, then let rank 0 handle the external GPU
  // calls, and finally broadcast the output vector.
  //
  // TODO: Ideally, we'd want to handle general architectures with multiple nodes and arbitrary
  // number of GPUs on each node. The parallel communications below are written using general vector
  // pull/push, so accounting for general architectures in the future should be trivial.

  // First, gather input vectors onto rank 0
  _console << "Serializing material " << getParam<MaterialName>("material") << std::endl;

  // Populate dofs to query
  std::unordered_map<processor_id_type, std::vector<dof_id_type>> query_dofs;
  if (processor_id() == 0)
    for (const auto & elem : _mesh.getMesh().element_ptr_range())
      if (elem->processor_id() != 0)
        for (auto qp : make_range(_material_data->nQPoints()))
          query_dofs[elem->processor_id()].push_back(dof(elem->id(), qp));

  // Gather the input dofs onto rank 0
  if (n_processors() > 1)
  {
    for (auto & it : _input)
      serialize(it.second, query_dofs);
    for (auto & it : _input_old)
      serialize(it.second, query_dofs);
  }

  // Now that input vectors are ready on rank 0, let's do the external calls to GPU(s) which are
  // supposed to be well-vectorized.
  _console << "Computing serialized material " << getParam<MaterialName>("material") << std::endl;
  if (processor_id() == 0)
    GPUCalls();

  // Now, the output vector is ready on rank 0, let's broadcast it to other processors
  _console << "Broadcasting computed material " << getParam<MaterialName>("material") << std::endl;
  _communicator.broadcast(query_dofs);
  if (n_processors() > 1)
    for (auto & it : _output)
      deserialize(it.second, query_dofs);

  _ready = true;
}

void
VectorizedMaterialBase::threadJoin(const UserObject & /*y*/)
{
}

void
VectorizedMaterialBase::allocateVector(std::map<std::string, std::vector<Real>> & container,
                                       const std::set<std::string> & prop_names)
{
  for (const auto & prop_name : prop_names)
  {
    if (_type[prop_name] == typeid(Real).name())
      allocateVectorHelper<Real>(container, prop_name);
    else if (_type[prop_name] == typeid(RankTwoTensor).name())
      allocateVectorHelper<RankTwoTensor>(container, prop_name);
    else if (_type[prop_name] == typeid(RankFourTensor).name())
      allocateVectorHelper<RankFourTensor>(container, prop_name);
    else
      mooseError("Unknown material property type for ", prop_name);
  }
}

void
VectorizedMaterialBase::fillVector(std::map<std::string, std::vector<Real>> & container,
                                   const std::set<std::string> & prop_names,
                                   const MaterialProperties & props)
{
  for (const auto & name : prop_names)
  {
    if (_type.at(name) == typeid(Real).name())
    {
      const auto & prop = getProp<Real>(props, name);
      for (auto qp : make_range(_material_data->nQPoints()))
        fill<Real>(container, name, _current_elem->id(), qp, prop[qp]);
    }
    if (_type.at(name) == typeid(RankTwoTensor).name())
    {
      const auto & prop = getProp<RankTwoTensor>(props, name);
      for (auto qp : make_range(_material_data->nQPoints()))
        fill<RankTwoTensor>(container, name, _current_elem->id(), qp, prop[qp]);
    }
    if (_type.at(name) == typeid(RankFourTensor).name())
    {
      const auto & prop = getProp<RankFourTensor>(props, name);
      for (auto qp : make_range(_material_data->nQPoints()))
        fill<RankFourTensor>(container, name, _current_elem->id(), qp, prop[qp]);
    }
  }
}

void
VectorizedMaterialBase::serialize(
    std::vector<Real> & vec,
    const std::unordered_map<processor_id_type, std::vector<dof_id_type>> & query_dofs)
{
  // Answer queries received from rank 0
  auto gather_data = [this, &vec](const processor_id_type /*pid*/,
                                  const std::vector<dof_id_type> & dofs,
                                  std::vector<Real> & data_to_fill)
  {
    data_to_fill.resize(dofs.size());
    for (auto i : index_range(dofs))
      data_to_fill[i] = vec[dofs[i]];
  };

  // Gather answers received from other processors
  auto act_on_data = [this, &vec](const processor_id_type /*pid*/,
                                  const std::vector<dof_id_type> & dofs,
                                  const std::vector<Real> & filled_data)
  {
    for (auto i : index_range(dofs))
      vec[dofs[i]] = filled_data[i];
  };

  const Real * datum = nullptr;
  libMesh::Parallel::pull_parallel_vector_data(
      _communicator, query_dofs, gather_data, act_on_data, datum);
}

void
VectorizedMaterialBase::deserialize(
    std::vector<Real> & vec,
    const std::unordered_map<processor_id_type, std::vector<dof_id_type>> & query_dofs)
{
  std::unordered_map<processor_id_type, std::vector<Real>> push_vecs;
  // Prepare data to push
  if (processor_id() == 0)
    for (const auto & [pid, pdofs] : query_dofs)
      for (auto pdof : pdofs)
        push_vecs[pid].push_back(vec[pdof]);

  // Act on data received from rank 0
  auto act_on_data =
      [this, &vec, &query_dofs](const processor_id_type /*pid*/, const std::vector<Real> & vec_recv)
  {
    const auto & dofs = query_dofs.at(processor_id());
    for (auto i : index_range(dofs))
      vec[dofs[i]] = vec_recv[i];
  };

  libMesh::Parallel::push_parallel_vector_data(_communicator, push_vecs, act_on_data);
}
