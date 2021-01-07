[![Build Status](https://github.com/salilab/optrep/workflows/build/badge.svg?branch=main)](https://github.com/salilab/optrep/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/salilab/optrep/branch/main/graph/badge.svg)](https://codecov.io/gh/salilab/pmi)

This repository is meant to be used as an [IMP](https://integrativemodeling.org) (*Integrative Modeling Platform*) module. 
It implements a method for obtaining an optimal coarse-grained bead representation for a given system, given input data, scoring functions, and sampling scheme, using an incremental coarse-graining algorithm (see publications below). 

## Note
See the wiki pages in this repo for details on how to use this module.\
Also see https://salilab.org/optimal_representation for examples using this module. 

## Prerequisites 
Note that a patch was applied to the `dof/__init__.py` file in IMP's PMI module (`pmi/pyext/src/`). The patched file can be found in [prereqs/pmi/pyext/src/dof/__init__.py](prereqs/pmi/pyext/src/dof/__init__.py). 

*Reason for the patch*: Monte Carlo move sizes vary based on bead size. Since we are now using non-uniform resolution beads instead of the earlier uniform resolution (fixed size) beads, we modify the Monte Carlo maximum translation parameter for each bead to be based on its size. Accordingly, we take as input a list of maximum translation parameters (`fbmaxtrans`), one for each coarse-grained bead size in the algorithm. For example, if we are coarse-graining incrementally to 1,5, and 10 residues per bead, we provide the maximum translation parameters for beads of 1,5, and 10 residues. For all other intervening bead sizes, we use cubic spline interpolation to find their maximum translation parameters. 

This patch will be committed to the master branch of IMP in a subsequent version, along with this module. 

## Information
_Author(s)_: Shruthi Viswanath 

_Maintainer_: `benmwebb`

_Date_: October 31st, 2018 

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Testable_: Yes.

_Parallelizeable_: Yes

_Publications_:
- S. Viswanath and A. Sali, Optimizing model representation for integrative structure determination of macromolecular assemblies, submitted. 
