# XFEMFunctionNeumannBC

!syntax description /DiracKernels/XFEMFunctionNeumannBC

## Overview

`XFEMFunctionNeumannBC` applies a Neumann BC at an interface cut by XFEM. The value of the BC is provided by a `Function`.

## Example Input Syntax

!listing test/tests/neumann_bc/2d_function.i block=DiracKernels

!syntax parameters /DiracKernels/XFEMFunctionNeumannBC

!syntax inputs /DiracKernels/XFEMFunctionNeumannBC

!syntax children /DiracKernels/XFEMFunctionNeumannBC
