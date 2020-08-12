/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                       BlackBear                              */
/*                                                              */
/*           (c) 2017 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#pragma once

#include "ADLAROMANCEStressUpdateBase.h"

class P91LAROMANCEStressUpdate : public ADLAROMANCEStressUpdateBase
{
public:
  static InputParameters validParams();

  P91LAROMANCEStressUpdate(const InputParameters & parameters);

protected:
  virtual std::vector<std::vector<std::vector<ROMInputTransform>>> getTransform() override;
  virtual std::vector<std::vector<std::vector<Real>>> getTransformCoefs() override;
  virtual std::vector<std::vector<std::vector<Real>>> getNormalizationLimits() override;
  virtual std::vector<std::vector<std::vector<Real>>> getInputLimits() override;
  virtual std::vector<std::vector<std::vector<Real>>> getCoefs() override;
  virtual std::vector<unsigned int> getTilings() override;
};
