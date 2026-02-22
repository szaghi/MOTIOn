---
title: About MOTIOn
---

# About MOTIOn

**MOTIOn** (Modular HPC Optimized Toolkit for IO in Fortran) is a pure Fortran 2003+ library for handling Input/Output operations in HPC scenarios for modern Fortran projects.

The main goal is to let the user deal with only a reduced set of operations while MOTIOn deals, in background, with the complexity of generating efficient parallel IO — exploiting HPC libraries like HDF5 with XDMF XML descriptors. The primary focus is on CFD solvers where a simulation domain is discretized into numerical grids and governing laws are integrated in parallel by MPI processes: MOTIOn provides a convenient API to IO those grids and fields with very low user effort.

## Authors

- Stefano Zaghi — [@szaghi](https://github.com/szaghi)

Contributions are welcome — see the [Contributing](contributing) page.

## Copyrights

MOTIOn is distributed under a multi-licensing system:

| Use case | License |
|----------|---------|
| FOSS projects | [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html) |
| Closed source / commercial | [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause) |
| Closed source / commercial | [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause) |
| Closed source / commercial | [MIT](http://opensource.org/licenses/MIT) |

> Anyone interested in using, developing, or contributing to MOTIOn is welcome — pick the license that best fits your needs.
