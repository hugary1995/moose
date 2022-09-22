//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ExternalWork.h"

registerMooseObject("TensorMechanicsApp", ExternalWork);

InputParameters
ExternalWork::validParams()
{
  InputParameters params = NodalPostprocessor::validParams();
  params.addClassDescription("This class computes the total external work. The power expenditure "
                             "(rate of external work) is defined as $\\mathcal{P}^\\text{ext} = "
                             "\\int_\\bodyboundary \\bft \\cdot \\dot{\\bs{\\phi}} \\diff{A}$. The "
                             "power expenditure is integrated in time to get the total work.");
  params.addRequiredCoupledVar("forces",
                               "The reaction forces associated with each of the displacement");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  return params;
}

ExternalWork::ExternalWork(const InputParameters & parameters)
  : NodalPostprocessor(parameters),
    _sum(0),
    _ndisp(coupledComponents("displacements")),
    _u(coupledValues("displacements")),
    _u_old(coupledValuesOld("displacements")),
    _nforce(coupledComponents("forces")),
    _forces(coupledValues("forces")),
    _forces_old(coupledValuesOld("forces")),
    _sum_old(getPostprocessorValueOldByName(name()))
{
  _u.resize(3, &_zero);
  _u_old.resize(3, &_zero);
  _forces.resize(3, &_zero);
  _forces_old.resize(3, &_zero);
}

void
ExternalWork::initialize()
{
  _sum = 0;
}

void
ExternalWork::execute()
{
  _sum += computeQpValue();
}

Real
ExternalWork::computeQpValue()
{
  RealVectorValue u((*_u[0])[_qp], (*_u[1])[_qp], (*_u[2])[_qp]);
  RealVectorValue u_old((*_u_old[0])[_qp], (*_u_old[1])[_qp], (*_u_old[2])[_qp]);
  RealVectorValue f((*_forces[0])[_qp], (*_forces[1])[_qp], (*_forces[2])[_qp]);
  RealVectorValue f_old((*_forces_old[0])[_qp], (*_forces_old[1])[_qp], (*_forces_old[2])[_qp]);

  return (f + f_old) / 2 * (u - u_old);
}

Real
ExternalWork::getValue()
{
  return _sum + _sum_old;
}

void
ExternalWork::finalize()
{
  gatherSum(_sum);
}

void
ExternalWork::threadJoin(const UserObject & y)
{
  const ExternalWork & pps = static_cast<const ExternalWork &>(y);
  _sum += pps._sum;
}
