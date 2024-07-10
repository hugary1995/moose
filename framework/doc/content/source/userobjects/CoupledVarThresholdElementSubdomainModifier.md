# CoupledVarThresholdElementSubdomainModifier

!syntax description /UserObjects/CoupledVarThresholdElementSubdomainModifier

## Overview

The `CoupledVarThresholdElementSubdomainModifier` can model

- Element death (with applications in ablation, fracture, etc.);
- Element activation (with applications in additive manufacturing, sintering, solidification, etc.);
- Moving interface (with applications in metal oxidation, phase transformation, melt pool, etc.).

The `CoupledVarThresholdElementSubdomainModifier` changes the element subdomain based on the given criterion. It also handles the corresponding

- Moving boundary/interface nodeset/sideset modification,
- Solution initialization, and
- Stateful material property initialization,

all of which are demonstrated using the following example.

Consider a unit square domain, and an auxiliary variable defined by the function $\phi(x,y,t) = (x-t)^2+y^2-0.5^2$. The function represents a signed-distance function of a circle of radius $0.5$ whose center is moving along the x-axis towards the right.

Initially, the domain is decomposed by a vertical line $x=0.25$. The elements on the left side of the vertical line have subdomain ID of 1, and the elements on the right side have subdomain ID of 2. The `CoupledVarThresholdElementSubdomainModifier` is used to change the subdomain ID from 2 to 1 for elements within the circle.

## Reversible vs. irreversible modification

If the `CoupledVarThresholdElementSubdomainModifier` is applied onto the entire domain, and the parameter `complement_subdomain_id` is set to 2, then the subdomain ID of all elements outside the circle will be set to 2:

!listing test/tests/userobjects/element_subdomain_modifier/reversible.i start=[moving_circle] end=[] include-end=true

!media media/userobjects/esm_reversible.jpg style=float:center;width:100%; caption=The result of a reversible element subdomain modifier at three different time steps

However, in many applications, e.g. element death and activation, the equivalent movement of element subdomains is not reversible. In this case, omitting the parameter `complement_subdomain_id` will make the subdomain modification irreversible:

!listing test/tests/userobjects/element_subdomain_modifier/irreversible.i start=[moving_circle] end=[] include-end=true

!media media/userobjects/esm_irreversible.jpg style=float:center;width:100%; caption=The result of an irreversible element subdomain modifier at three different time steps

## Moving boundary/interface nodeset/sideset modification

The change of element subdomains will alter the definition of certain sidesets and nodesets. The `CoupledVarThresholdElementSubdomainModifier` optionally takes the parameter `moving_boundaries` to help modify the corresponding sideset/nodeset. If the boundary provided through the `moving_boundaries` parameter already exists, the modifier will attempt to modify the provided sideset/nodeset whenever an element changes subdomain. If the boundary does not exist, the modifier will create a sideset and a nodeset with the provided name. This parameter must also be set with `moving_boundary_subdomain_pairs`, to specify the subdomains the boundary lies between. 

!listing test/tests/userobjects/element_subdomain_modifier/moving_boundary.i start=[moving_circle] end=[] include-end=true

!media media/userobjects/esm_sideset.jpg style=float:center;width:100%; caption=The evolving sideset (green) at three different time steps, without a provided boundary from the user.

!media media/userobjects/esm_nodeset.jpg style=float:center;width:100%; caption=The evolving nodeset (green) at three different time steps, without a provided boundary from the user.

The parameter `moving_boundary_subdomain_pairs` also allows for the external parts of the geometry to be included in the moving boundary by setting the subdomain pair to only one subdomain.

!listing test/tests/userobjects/element_subdomain_modifier/external_moving_boundary.i start=[moving_circle] end=[] include-end=true

Nodal and integrated BCs can be applied on the moving boundary.

## Solution initialization

Depending on the physics, one may or may not want to reinitialize the solution when an element and its related nodes change subdomain. The parameter `subdomains_to_reinitialize` defaults to the whole geometry, thus reinitializing any changed element. 

!listing test/tests/userobjects/element_subdomain_modifier/reinitialization.i start=[UserObjects] end=[AuxVariables]

If a list of subdomains is set for this parameter, then only elements which change into a subdomain in the given list will be reinitialized. If an empty list is given, then no changed elements will be reinitialized.

!listing test/tests/userobjects/element_subdomain_modifier/reinitialization_into.i start=[UserObjects] end=[AuxVariables]

A further restriction can be implemented, where only elements which begin outside of `subdomains_to_reinitialize` and then change into those subdomains, will be reinitialized. This has applications for element activation, where only elements moving from inactive to active domains need to be reinitialized.

!listing test/tests/userobjects/element_subdomain_modifier/reinitialization_from_into.i start=[UserObjects] end=[AuxVariables]

Suppose initially there is an auxiliary variable $u=1$ everywhere inside the domain, and the variable value in subdomain 1 (blue) doubles at each time step:

!listing test/tests/userobjects/element_subdomain_modifier/initial_condition.i start=[AuxVariables] end=[Postprocessors]

!media media/userobjects/esm_ic.jpg style=float:center;width:100%; caption=The auxiliary variable $u$ at three different time steps

## Stateful material property initialization

Similarly, all stateful material properties will be re-initialized when an element changes subdomain. Suppose initially the diffusivity is $0.5$ everywhere, and the diffusivity doubles at each time step:

!listing test/tests/userobjects/element_subdomain_modifier/stateful_property.i start=[Materials] end=[Executioner]

!media media/userobjects/esm_material.jpg style=float:center;width:100%; caption=The diffusivity at three different time steps

!syntax parameters /UserObjects/CoupledVarThresholdElementSubdomainModifier

!syntax inputs /UserObjects/CoupledVarThresholdElementSubdomainModifier

!syntax children /UserObjects/CoupledVarThresholdElementSubdomainModifier
