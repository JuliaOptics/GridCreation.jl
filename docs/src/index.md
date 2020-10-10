# GridCreation

GridCreation includes simple functions for producing 2D grids on which to compute functions.  Cartesian grids are only produced as vectors, as broadcasting allows them to be used as a grid with lower memory allocation.  Polar grids are generated as dense 2D arrays due to their lack of exploitable symmetry in a cartesian frame of reference.

## Installation

GridCreation is a registered Julia package, it may be installed the usual way:
```sh
julia> Pkg.add("GridCreation")
```

## Rationale

In Julia, it is viewed as idiomatic to use "meshgrid" for Cartesian grids.  This is because dot syntax broadcasts "correctly" between row and column vectors, so the grid need not be stored.  For this reason, mkCartVecs produces only vectors.  Because these is no "mixing" between X and Y, there is not a more performant approach.

cartVecsToPolarGrid explicitly creates a meshgrid for its output.  This is because both the radial and azimuthal variables require some computation.  The duplicate computation is more expensive than the allocation when the grid is used several times, as in optical simulation and analysis.  For this reason, cartVecsToPolarGrid forms a grid.

## Functions
```@docs
mkCartVecs
cartVecsToPolarGrid
```
