# GridCreation

GridCreation includes simple functions for producing 2D grids on which to compute functions.  Cartesian grids are only produced as vectors, as broadcasting allows them to be used as a grid with lower memory allocation.  Polar grids are generated as dense 2D arrays due to their lack of exploitable symmetry in a cartesian frame of reference.

## Installation

GridCreation is a registered Julia package, it may be installed the usual way:
```sh
julia> Pkg.add("GridCreation")
```

## Functions
```@docs
mkCartVecs
cartVecsToPolarGrid
```
