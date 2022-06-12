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

#include <any>

/* This class gathers material properties requested by a material object to form a vector */
class VectorizedMaterial : public ElementUserObject, public DependencyResolverInterface
{
public:
  static InputParameters validParams();

  VectorizedMaterial(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject &) override;
  virtual const std::set<std::string> & getRequestedItems() override { return _reqs; }
  virtual const std::set<std::string> & getSuppliedItems() override { return _sups; }

  virtual Real getReal(const std::string name, unsigned int elem, unsigned int qp) const;
  virtual RankTwoTensor
  getRankTwoTensor(const std::string name, unsigned int elem, unsigned int qp) const;
  virtual RankFourTensor
  getRankFourTensor(const std::string name, unsigned int elem, unsigned int qp) const;

  virtual bool ready() const { return _ready; }

protected:
  MaterialBase * _mat;
  const std::vector<MaterialPropertyName> _mat_prop_names;
  const std::vector<std::string> _mat_prop_types;
  std::map<MaterialPropertyName, std::string> _type;

  std::map<std::string, unsigned int> _prop_name_id_map;
  std::map<unsigned int, std::string> _prop_id_name_map;

  mutable std::set<std::string> _reqs;
  std::set<std::string> _sups;

  std::set<std::string> _input;
  std::set<std::string> _input_old;
  std::set<std::string> _output;

  std::map<std::string, std::vector<MooseArray<Real>>> _input_Real;
  std::map<std::string, std::vector<MooseArray<RankTwoTensor>>> _input_RankTwoTensor;
  std::map<std::string, std::vector<MooseArray<RankFourTensor>>> _input_RankFourTensor;

  std::map<std::string, std::vector<MooseArray<Real>>> _input_Real_old;
  std::map<std::string, std::vector<MooseArray<RankTwoTensor>>> _input_RankTwoTensor_old;
  std::map<std::string, std::vector<MooseArray<RankFourTensor>>> _input_RankFourTensor_old;

  std::map<std::string, std::vector<MooseArray<Real>>> _output_Real;
  std::map<std::string, std::vector<MooseArray<RankTwoTensor>>> _output_RankTwoTensor;
  std::map<std::string, std::vector<MooseArray<RankFourTensor>>> _output_RankFourTensor;

private:
  void allocateInputVector(const std::string & name);
  void allocateOutputVector(const std::string & name);
  void allocateInputOldVector(const std::string & name);
  bool _ready;
};
