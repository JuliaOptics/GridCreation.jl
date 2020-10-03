module GridCreation

export CFFT, CInterSample, mkCartVecs, cartVecsToPolarGrid

const CFFT = -1;
const CInterSample = -2;

# this function is lifted from pkg Fourier.
# A little copying is better than a little dependency.
function _fftrange(n)
    if iseven(n)
        return -(n ÷ 2):(n÷2-1)
    end
    return -(n ÷ 2):(n÷2)
end
"""
    mkCartVecs(dx, size[, center])

Create a 2D "grid" by returning a row vector x and column vector y.

The grid will have sample spacing dx, size[0] == size(x), and size[1] == size(y).

If size is an int, it is broadcast to [size,size].

The grid will be centered according to such that the center element contains zero.  Special
values for center include -1 (default, FFT-like center) and -2, for "interpixel"
sampling, where 0 is at size/2.  Otherwise, zero is at the centerth element.

The centering rule may not be different between x and y.

# Examples

## FFT-aligned grid

Note that the ceil rounded center element contains zero.
```julia-repl
julia> x,y = mkCartVecs(.1, 4); # implicit center=CFFT
julia> x
1×4 LinearAlgebra.Adjoint{Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:
 -0.2  -0.1  0.0  0.1
julia> collect(y)
4-element Array{Float64,1}:
 -0.2
 -0.1
  0.0
  0.1
```

## Intersample aligned grid
```julia-repl
x,y=mkCartVecs(.1, 4, center=CInterSample);

julia> x
1×4 LinearAlgebra.Adjoint{Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:
 -0.15  -0.05  0.05  0.15
```

## index-centered grid
```
julia> x,y=mkCartVecs(.1, 4, center=1);
julia> x
1×4 LinearAlgebra.Adjoint{Float64,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}:
 1.11022e-17  0.1  0.2  0.3
```
"""
function mkCartVecs(dx, size::Tuple{Integer,Integer}; center::Integer=CFFT)
    # fftrange produces a range object.  dot add shifts the range
    # and multiplication scales it.
    if center == CFFT
        X = (_fftrange(size[2])*dx)';
        Y = _fftrange(size[1])*dx;
    elseif center == CInterSample
        X = (_fftrange(size[2])*dx.+(dx/2))';
        Y = (_fftrange(size[1])*dx.+(dx/2));
    else
        X = _fftrange(size[2])*dx;
        Y = _fftrange(size[1])*dx;
        X = (X .- X[center])';
        Y = (Y .- Y[center]);
    end
    return X, Y;
end

function mkCartVecs(dx, size::Integer; center::Integer=CFFT)
    return mkCartVecs(dx, (size,size), center=center);
end

"""
    cartVecsToPolarGrid(X, Y)

Construct a polar grid (ρ,θ) from the row vector X and column vector Y.  Returns
a pair of 2D arrays containing the radial and azimuthal coordinates.

# Examples
```julia-repl
julia> X=(-2:1)'; Y=-2:1;
julia> ρ,θ=cartVecsToPolarGrid(X,Y);
julia> ρ
4×4 Array{Float64,2}:
 2.82843  2.23607  2.0  2.23607
 2.23607  1.41421  1.0  1.41421
 2.0      1.0      0.0  1.0
 2.23607  1.41421  1.0  1.41421
julia> θ
4×4 Array{Float64,2}:
 -2.35619  -2.03444  -1.5708  -1.10715
 -2.67795  -2.35619  -1.5708  -0.785398
  3.14159   3.14159   0.0      0.0
  2.67795   2.35619   1.5708   0.785398
```
See also: [`mkCartVecs`](@ref)
"""
function cartVecsToPolarGrid(X, Y)
    # I don't think allocating a meshgrid is avoidable here.
    # the rho matrix can be made with a dot on x and y,
    # but we would need a 2D view into both X and Y for atan
    # grid_x = [x for x in X, y in Y;
    # grid_y = [y for y in Y, x in X;
    ρ = @. sqrt(X^2 + Y^2);
    θ = [atan(y,x) for y in Y, x in X'];
    return ρ, θ
end

end
