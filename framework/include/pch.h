//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// metaphysicl
#include "metaphysicl/dualnumber_decl.h"
#include "metaphysicl/ct_types.h"
#include "metaphysicl/raw_type.h"
#include "metaphysicl/dualsemidynamicsparsenumberarray.h"
#include "metaphysicl/dynamic_std_array_wrapper.h"
#include "metaphysicl/dualnumberarray.h"

// libmesh
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh.h"
#include "libmesh/remote_elem.h"
#include "libmesh/elem.h"
#include "libmesh/elem_range.h"
#include "libmesh/mesh_base.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/bounding_box.h"
#include "libmesh/point.h"
#include "libmesh/type_vector.h"
#include "libmesh/tensor_tools.h"

// contrib
#include "hit/hit.h"
#include "hit/parse.h"
#include "nlohmann/json.h"

// actions
#include "Action.h"

// auxkernels
#include "AuxKernel.h"
#include "AuxKernelBase.h"

// auxscalarkernels
#include "AuxScalarKernel.h"

// base
#include "MooseObject.h"
#include "ResidualObject.h"
#include "Moose.h"
#include "libMeshReducedNamespace.h"
#include "Assembly.h"
#include "MooseApp.h"
#include "ParallelParamObject.h"
#include "MooseBase.h"
#include "Adaptivity.h"
#include "MooseError.h"
#include "NeighborResidualObject.h"
#include "Registry.h"
#include "Factory.h"
#include "MooseFunctor.h"
#include "TheWarehouse.h"

// bcs
#include "BoundaryCondition.h"
#include "IntegratedBCBase.h"
#include "IntegratedBC.h"
#include "ADIntegratedBC.h"
#include "NodalBCBase.h"
#include "NodalBC.h"
#include "ADNodalBC.h"

// constraints
#include "Constraint.h"

// controls
#include "Control.h"
#include "ControllableItem.h"
#include "ControllableParameter.h"

// distributions
#include "DistributionInterface.h"

// functions
#include "Function.h"
#include "FunctionInterface.h"

// functormaterials
#include "FunctorMaterialProperty.h"

// fvkernels
#include "FVKernel.h"

// ics
#include "InitialConditionWarehouse.h"
#include "InitialConditionTempl.h"

// interfaces
#include "TaggingInterface.h"
#include "Coupleable.h"
#include "ElementIDInterface.h"
#include "MeshMetaDataInterface.h"
#include "DependencyResolverInterface.h"
#include "FunctorInterface.h"
#include "BlockRestrictable.h"
#include "BoundaryRestrictable.h"
#include "MeshChangedInterface.h"
#include "OutputInterface.h"
#include "SetupInterface.h"
#include "TransientInterface.h"

// kernels
#include "KernelBase.h"
#include "Kernel.h"
#include "ADKernel.h"

// materials
#include "Material.h"
#include "MaterialBase.h"
#include "DerivativeMaterialInterface.h"
#include "MaterialData.h"
#include "MaterialProperty.h"
#include "MaterialPropertyInterface.h"

// mesh
#include "MooseMesh.h"

// meshgenerators
#include "MeshGenerator.h"

// parser
#include "Builder.h"
#include "Parser.h"
#include "Syntax.h"

// postprocessors
#include "GeneralPostprocessor.h"
#include "Postprocessor.h"
#include "ElementPostprocessor.h"

// problems
#include "FEProblemBase.h"
#include "SubProblem.h"
#include "FEProblem.h"

// reporters
#include "ReporterData.h"
#include "Reporter.h"
#include "GeneralReporter.h"
#include "ReporterInterface.h"
#include "ReporterContext.h"

// restart
#include "Restartable.h"
#include "RestartableData.h"
#include "DataIO.h"
#include "RestartableEquationSystems.h"

// samplers
#include "SamplerInterface.h"

// systems
#include "SystemBase.h"
#include "NonlinearSystemBase.h"
#include "AuxiliarySystem.h"
#include "NonlinearSystem.h"

// userobjects
#include "UserObject.h"
#include "GeneralUserObject.h"
#include "ElementUserObject.h"
#include "SideUserObject.h"
#include "ThreadedGeneralUserObject.h"

// utils
#include "MooseTypes.h"
#include "MooseUtils.h"
#include "InputParameters.h"
#include "MooseArray.h"
#include "EigenADReal.h"
#include "RankTwoTensor.h"
#include "RankTwoTensorImplementation.h"
#include "RankFourTensor.h"
#include "RankFourTensorImplementation.h"
#include "StreamArguments.h"
#include "ColumnMajorMatrix.h"
#include "MathUtils.h"

// variables
#include "MooseVariableFE.h"
#include "MooseVariableField.h"
#include "MooseVariableData.h"
#include "VariableWarehouse.h"

// warehouses
#include "MooseObjectWarehouse.h"
#include "ExecuteMooseObjectWarehouse.h"
