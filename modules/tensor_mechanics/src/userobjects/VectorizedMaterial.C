//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VectorizedMaterial.h"
#include "MooseMesh.h"

registerMooseObject("TensorMechanicsApp", VectorizedMaterial);

InputParameters
VectorizedMaterial::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addRequiredParam<MaterialName>("material", "The material to vectorize.");
  return params;
}

VectorizedMaterial::VectorizedMaterial(const InputParameters & parameters)
  : ElementUserObject(parameters), _ready(false)
{
  _sups.insert(name());
}

void
VectorizedMaterial::initialSetup()
{
  _mat = &getMaterialByName(getParam<MaterialName>("material"), /*no_warn=*/true);

  const auto & storage = _material_data->getMaterialPropertyStorage();
  const auto & prop_name_id_map = storage.propIDs();
  for (const auto & [name, id] : prop_name_id_map)
  {
    _prop_name_id_map.emplace(name, id);
    _prop_id_name_map.emplace(id, name);
    _type.emplace(name, _material_data->props()[id]->type());
  }

  const auto & reqs = _mat->getRequestedItems();
  const auto & sups = _mat->getSuppliedItems();
  const auto & dep_ids = _mat->getMatPropDependencies();
  _input.insert(reqs.begin(), reqs.end());
  _output.insert(sups.begin(), sups.end());
  for (const auto & dep_id : dep_ids)
    if (storage.isStatefulProp(_prop_id_name_map[dep_id]))
      _input_old.insert(_prop_id_name_map[dep_id]);

  // The requested material properties must be reinited before the execution of this userobject.
  _reqs.insert(reqs.begin(), reqs.end());
  _material_property_dependencies.insert(dep_ids.begin(), dep_ids.end());

  // Debugging output
  _console << "Input material property names: ";
  for (const auto & name : _input)
    _console << name << "[" << _type[name] << "] ";
  _console << std::endl;
  _console << "Input old material property names: ";
  for (const auto & name : _input_old)
    _console << name << "[" << _type[name] << "] ";
  _console << std::endl;
  _console << "Output material property names: ";
  for (const auto & name : _output)
    _console << name << "[" << _type[name] << "] ";
  _console << std::endl;

  for (const auto & name : _input)
    allocateInputVector(name);
  for (const auto & name : _output)
    allocateOutputVector(name);
  for (const auto & name : _input_old)
    allocateInputOldVector(name);
}

void
VectorizedMaterial::allocateInputVector(const std::string & name)
{
  try
  {
    unsigned int N = _mesh.maxElemId();
    if (_type.at(name) == typeid(Real).name())
      _input_Real[name] = std::vector<MooseArray<Real>>(N);
    if (_type.at(name) == typeid(RankTwoTensor).name())
      _input_RankTwoTensor[name] = std::vector<MooseArray<RankTwoTensor>>(N);
    if (_type.at(name) == typeid(RankFourTensor).name())
      _input_RankFourTensor[name] = std::vector<MooseArray<RankFourTensor>>(N);
  }
  catch (std::out_of_range &)
  {
    mooseError("Unknown material property type for ", name);
  }
}

void
VectorizedMaterial::allocateOutputVector(const std::string & name)
{
  try
  {
    unsigned int N = _mesh.maxElemId();
    if (_type.at(name) == typeid(Real).name())
      _output_Real[name] = std::vector<MooseArray<Real>>(N);
    if (_type.at(name) == typeid(RankTwoTensor).name())
      _output_RankTwoTensor[name] = std::vector<MooseArray<RankTwoTensor>>(N);
    if (_type.at(name) == typeid(RankFourTensor).name())
      _output_RankFourTensor[name] = std::vector<MooseArray<RankFourTensor>>(N);
  }
  catch (std::out_of_range &)
  {
    mooseError("Unknown material property type for ", name);
  }
}

void
VectorizedMaterial::allocateInputOldVector(const std::string & name)
{
  try
  {
    unsigned int N = _mesh.maxElemId();
    if (_type.at(name) == typeid(Real).name())
      _input_Real_old[name] = std::vector<MooseArray<Real>>(N);
    if (_type.at(name) == typeid(RankTwoTensor).name())
      _input_RankTwoTensor_old[name] = std::vector<MooseArray<RankTwoTensor>>(N);
    if (_type.at(name) == typeid(RankFourTensor).name())
      _input_RankFourTensor_old[name] = std::vector<MooseArray<RankFourTensor>>(N);
  }
  catch (std::out_of_range &)
  {
    mooseError("Unknown material property type for ", name);
  }
}

void
VectorizedMaterial::initialize()
{
  _ready = false;
}

void
VectorizedMaterial::execute()
{
  for (const auto & name : _input)
  {
    auto prop_id = _prop_name_id_map[name];
    if (_type.at(name) == typeid(Real).name())
      _input_Real[name][_current_elem->id()] =
          static_cast<MaterialPropertyBase<Real, false> *>(_material_data->props()[prop_id])->get();
    if (_type.at(name) == typeid(RankTwoTensor).name())
      _input_RankTwoTensor[name][_current_elem->id()] =
          static_cast<MaterialPropertyBase<RankTwoTensor, false> *>(
              _material_data->props()[prop_id])
              ->get();
    if (_type.at(name) == typeid(RankFourTensor).name())
      _input_RankFourTensor[name][_current_elem->id()] =
          static_cast<MaterialPropertyBase<RankFourTensor, false> *>(
              _material_data->props()[prop_id])
              ->get();
  }

  for (const auto & name : _input_old)
  {
    auto prop_id = _prop_name_id_map[name];
    if (_type.at(name) == typeid(Real).name())
      _input_Real_old[name][_current_elem->id()] =
          static_cast<MaterialPropertyBase<Real, false> *>(_material_data->props()[prop_id])->get();
    if (_type.at(name) == typeid(RankTwoTensor).name())
      _input_RankTwoTensor_old[name][_current_elem->id()] =
          static_cast<MaterialPropertyBase<RankTwoTensor, false> *>(
              _material_data->props()[prop_id])
              ->get();
    if (_type.at(name) == typeid(RankFourTensor).name())
      _input_RankFourTensor_old[name][_current_elem->id()] =
          static_cast<MaterialPropertyBase<RankFourTensor, false> *>(
              _material_data->props()[prop_id])
              ->get();
  }

  for (const auto & name : _output)
  {
    if (_type.at(name) == typeid(Real).name())
      _output_Real[name][_current_elem->id()].resize(_material_data->nQPoints());
    if (_type.at(name) == typeid(RankTwoTensor).name())
      _output_RankTwoTensor[name][_current_elem->id()].resize(_material_data->nQPoints());
    if (_type.at(name) == typeid(RankFourTensor).name())
      _output_RankFourTensor[name][_current_elem->id()].resize(_material_data->nQPoints());
  }
}

void
VectorizedMaterial::finalize()
{
  // Fake the behavior of an external vectorized library
  for (auto elem : make_range(_mesh.maxElemId()))
    for (auto qp : make_range(_material_data->nQPoints()))
    {
      const RankFourTensor C = _input_RankFourTensor["elasticity_tensor"][elem][qp];
      const RankTwoTensor E = _input_RankTwoTensor["mechanical_strain"][elem][qp];
      _output_RankTwoTensor["small_stress"][elem][qp] = C * E;
      _output_RankFourTensor["small_jacobian"][elem][qp] = C;
    }

  _ready = true;
  _console << "Vectorized material " << getParam<MaterialName>("material") << std::endl;
}

void
VectorizedMaterial::threadJoin(const UserObject & /*y*/)
{
}

Real
VectorizedMaterial::getReal(const std::string name, unsigned int elem, unsigned int qp) const
{
  return _output_Real.at(name)[elem][qp];
}

RankTwoTensor
VectorizedMaterial::getRankTwoTensor(const std::string name,
                                     unsigned int elem,
                                     unsigned int qp) const
{
  return _output_RankTwoTensor.at(name)[elem][qp];
}

RankFourTensor
VectorizedMaterial::getRankFourTensor(const std::string name,
                                      unsigned int elem,
                                      unsigned int qp) const
{
  return _output_RankFourTensor.at(name)[elem][qp];
}
