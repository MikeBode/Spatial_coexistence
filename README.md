# Spatial_coexistence
MATLAB code to support "Complex patch geometry can promote species coexistence through a competition-colonisation trade-off". 

This paper can be found on biorXiv:

Explaining the maintenance of diverse species assemblages is a central goal of ecology and conservation. Recent coexistence mechanisms highlight the role of dispersal as a source of the differences that allow similar species to coexist. 

We have proposed a new mechanism for species coexistence that is based on dispersal differences, and on the geometry of the habitat patch. In a finite habitat patch, species with different dispersal abilities will arrange themselves in stable, concentric patterns of dominance. Species with superior competitive and dispersal abilities will dominate the interior of the patch, with inferior species at the periphery.

The code in this repository demonstrates this mechanism. One piece of code simulates the coexistence of two species in a 1D domain. A second piece of code simulates coexistence in 2D habitat patches with realistic geometries. A final set of code uses metrics from landscape ecology to demonstrate that habitat patches with more complex geometries can more easily support coexistence.

The code to create Figure 1 is called "Make_Figure1.m". It can be run as a standalone function. It requires Matlab's elliptical PDE solver, PDEPE.

The code to create Figure 2 is called "Make_Figure2.m". It first requires that the user run the code "Simulate_Island_equilibria_Fig_2.m".

The code to create Figure 3 is called "Make_Figure3.m". It first requires that the user run the code "Measure_coexistence_potential.m" for the desired number of random communities.
