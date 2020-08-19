//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

// libMesh include files.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/exodusII_io_helper.h"

#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <ctime>
#include <random>

// Bring in everything from the libMesh namespace
using namespace libMesh;

std::vector<Real> sample_gaussian(const std::vector<Real> & eigvals,
                                  const std::vector<std::vector<Real>> & eigvecs,
                                  std::default_random_engine & generator);

void compute_field(std::vector<Real> & P, const std::vector<Real> & Xi, Real mean, Real CV);

int
main(int argc, char ** argv)
{
  // Check for proper usage.
  if (argc != 1)
    libmesh_error_msg("\nUsage: " << argv[0]);

  // Initialize libMesh and the dependent libraries.
  LibMeshInit init(argc, argv);

  if (init.comm().size() > 1)
    libmesh_error_msg("Parallel sampling is not supported, use 1 processor only.");

  // Read the mesh
  Mesh mesh(init.comm());
  mesh.read("basis.e");

  // Read the eigenpairs
  std::vector<Real> eigvals;
  std::vector<std::vector<Real>> eigvecs;

  ExodusII_IO basis(mesh);
  basis.read("basis.e");
  ExodusII_IO_Helper & basis_helper = basis.get_exio_helper();
  for (int i = 1; i < basis.get_num_time_steps(); i++)
  {
    // read eigenvalue
    std::vector<Real> eigval;
    basis.read_global_variable({"d"}, i, eigval);
    eigvals.push_back(eigval[0]);

    // read eigenvector
    basis_helper.read_nodal_var_values("v", i);
    eigvecs.push_back(basis_helper.nodal_var_values);
  }

  // sample Gaussian fields
  std::default_random_engine generator;
  generator.seed(std::time(NULL));
  std::vector<Real> Xi = sample_gaussian(eigvals, eigvecs, generator);

  // transform to marginal Gamma fields
  std::vector<Real> field;
  compute_field(field, Xi, 14, 0.1);

  // write random field
  ExodusII_IO out(mesh);
  out.write("field.e");
  ExodusII_IO_Helper & out_helper = out.get_exio_helper();
  out_helper.initialize_nodal_variables({"field"});
  out_helper.write_nodal_values(1, field, 1);

  return EXIT_SUCCESS;
}

std::vector<Real>
sample_gaussian(const std::vector<Real> & eigvals,
                const std::vector<std::vector<Real>> & eigvecs,
                std::default_random_engine & generator)
{
  unsigned int ndof = eigvecs[0].size();
  std::vector<Real> Xi(ndof);
  std::normal_distribution<Real> distribution(0.0, 1.0);

  for (unsigned int i = 0; i < eigvals.size(); i++)
  {
    Real eta = distribution(generator);
    for (unsigned int j = 0; j < ndof; j++)
      Xi[j] += std::sqrt(eigvals[i]) * eta * eigvecs[i][j];
  }

  return Xi;
}

void
compute_field(std::vector<Real> & P, const std::vector<Real> & Xi, Real mean, Real CV)
{
  unsigned int ndof = Xi.size();

  // Normal distribution
  auto normal = boost::math::normal_distribution<Real>(0, 1);

  // Target distribution
  Real sigma = CV * mean;

  // Transform into the first field
  P.resize(ndof);
  for (unsigned int i = 0; i < ndof; i++)
  {
    Real p = 2 * boost::math::cdf(normal, Xi[i]);
    std::cout << p << ", ";
    if (p < 0 or p > 2)
      std::cout << "wrong wrong" << std::endl;
    P[i] = mean -
           sigma * std::sqrt(2) * boost::math::erfc_inv<Real>(2 * boost::math::cdf(normal, Xi[i]));
    std::cout << P[i] << std::endl;
  }
}
