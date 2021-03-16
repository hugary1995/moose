# XFEMConvectiveHeatFluxBC

!syntax description /DiracKernels/XFEMConvectiveHeatFluxBC

## Overview

`XFEMConvectiveHeatFluxBC` applies convective heat flux at an interface cut by XFEM. This is implemented with DiracKernel. This boundary condition computes convective heat flux $q'' = H \cdot (T - T_{inf})$, where $H$ is convective heat transfer coefficient,
$T$ is the temperature, and $T_{inf}$ is far field temperature.  Both $H$ and $T_{inf}$ are coupled as material properties.

## Example Input Syntax

!listing test/tests/heat_flux_bc/2d.i block=DiracKernels

!syntax parameters /DiracKernels/XFEMConvectiveHeatFluxBC

!syntax inputs /DiracKernels/XFEMConvectiveHeatFluxBC

!syntax children /DiracKernels/XFEMConvectiveHeatFluxBC
