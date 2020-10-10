var documenterSearchIndex = {"docs":
[{"location":"#GridCreation","page":"GridCreation","title":"GridCreation","text":"","category":"section"},{"location":"","page":"GridCreation","title":"GridCreation","text":"GridCreation includes simple functions for producing 2D grids on which to compute functions.  Cartesian grids are only produced as vectors, as broadcasting allows them to be used as a grid with lower memory allocation.  Polar grids are generated as dense 2D arrays due to their lack of exploitable symmetry in a cartesian frame of reference.","category":"page"},{"location":"#Installation","page":"GridCreation","title":"Installation","text":"","category":"section"},{"location":"","page":"GridCreation","title":"GridCreation","text":"GridCreation is a registered Julia package, it may be installed the usual way:","category":"page"},{"location":"","page":"GridCreation","title":"GridCreation","text":"julia> Pkg.add(\"GridCreation\")","category":"page"},{"location":"#Rationale","page":"GridCreation","title":"Rationale","text":"","category":"section"},{"location":"","page":"GridCreation","title":"GridCreation","text":"In Julia, it is viewed as idiomatic to use \"meshgrid\" for Cartesian grids.  This is because dot syntax broadcasts \"correctly\" between row and column vectors, so the grid need not be stored.  For this reason, mkCartVecs produces only vectors.  Because these is no \"mixing\" between X and Y, there is not a more performant approach.","category":"page"},{"location":"","page":"GridCreation","title":"GridCreation","text":"cartVecsToPolarGrid explicitly creates a meshgrid for its output.  This is because both the radial and azimuthal variables require some computation.  The duplicate computation is more expensive than the allocation when the grid is used several times, as in optical simulation and analysis.  For this reason, cartVecsToPolarGrid forms a grid.","category":"page"},{"location":"#Functions","page":"GridCreation","title":"Functions","text":"","category":"section"},{"location":"","page":"GridCreation","title":"GridCreation","text":"mkCartVecs\ncartVecsToPolarGrid","category":"page"},{"location":"#GridCreation.mkCartVecs","page":"GridCreation","title":"GridCreation.mkCartVecs","text":"mkCartVecs(dx, size[, center])\n\nCreate a 2D \"grid\" by returning a row vector x and column vector y.\n\nThe grid will have sample spacing dx, size[0] == size(x), and size[1] == size(y).\n\nIf size is an int, it is broadcast to [size,size].\n\nThe grid will be centered according to such that the center element contains zero.  Special values for center include -1 (default, FFT-like center) and -2, for \"interpixel\" sampling, where 0 is at size/2.  Otherwise, zero is at the centerth element.\n\nThe centering rule may not be different between x and y.\n\nExamples\n\nFFT-aligned grid\n\nNote that the ceil rounded center element contains zero.\n\njulia> x,y = mkCartVecs(.1, 4); # implicit center=CFFT\njulia> x\n1×4 LinearAlgebra.Adjoint{Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:\n -0.2  -0.1  0.0  0.1\njulia> collect(y)\n4-element Array{Float64,1}:\n -0.2\n -0.1\n  0.0\n  0.1\n\nIntersample aligned grid\n\nx,y=mkCartVecs(.1, 4, center=CInterSample);\n\njulia> x\n1×4 LinearAlgebra.Adjoint{Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:\n -0.15  -0.05  0.05  0.15\n\nindex-centered grid\n\njulia> x,y=mkCartVecs(.1, 4, center=1);\njulia> x\n1×4 LinearAlgebra.Adjoint{Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:\n 1.11022e-17  0.1  0.2  0.3\n\n\n\n\n\n","category":"function"},{"location":"#GridCreation.cartVecsToPolarGrid","page":"GridCreation","title":"GridCreation.cartVecsToPolarGrid","text":"cartVecsToPolarGrid(X, Y)\n\nConstruct a polar grid (ρ,θ) from the row vector X and column vector Y.  Returns a pair of 2D arrays containing the radial and azimuthal coordinates.\n\nExamples\n\njulia> X=(-2:1)'; Y=-2:1;\njulia> ρ,θ=cartVecsToPolarGrid(X,Y);\njulia> ρ\n4×4 Array{Float64,2}:\n 2.82843  2.23607  2.0  2.23607\n 2.23607  1.41421  1.0  1.41421\n 2.0      1.0      0.0  1.0\n 2.23607  1.41421  1.0  1.41421\njulia> θ\n4×4 Array{Float64,2}:\n -2.35619  -2.03444  -1.5708  -1.10715\n -2.67795  -2.35619  -1.5708  -0.785398\n  3.14159   3.14159   0.0      0.0\n  2.67795   2.35619   1.5708   0.785398\n\nSee also: mkCartVecs\n\n\n\n\n\n","category":"function"}]
}