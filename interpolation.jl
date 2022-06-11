function interpolate1D(f::Vector, x::Vector, x0)
    """Returns value of function f at x0 by linear interpolation of provided datapoints.
    x should be ascending. Interpolation is periodic"""

    bounds = getBoundsIndices(x, x0)
    i = bounds[1]
    j = bounds[2]

    # Interpolate and return
    f[i] + (f[j] - f[i])*(x0 - x[i])/(x[j] - x[i])
end

function interpolate2D(f::Matrix, x::Vector, y::Vector, x0, y0)
    """
    Returns value of function f at (x0,y0) by repeated linear interpolation of provided datapoints.
    x should ascend as row index increases. y should ascend as column index increases.
    Interpolation is periodic.
    """

    # Find upper and lower y indices of grid square containing (x0, y0)
    # We don't need the x indices
    bounds = getBoundsIndices(y, y0)
    i = bounds[1]
    j = bounds[2]

    # Intepolate for x at borders of grid (fixed y)
    fi = interpolate1D(f[:,i], x, x0)
    fj = interpolate1D(f[:,j], x, x0)

    return interpolate1D([fi, fj], [y[i], y[j]], y0)
end

function interpolate3D(f::Array{T, 3}, x::Vector, y::Vector, z::Vector, x0, y0, z0) where T
    """Returns value of function f at (x0,y0,z0) by bilinear interpolation of provided datapoints.
    x should be ascending. Interpolation is periodic."""

    # Find bounding incides of grid cube containing (x0, y0, z0)\
    # We don't need the x indices
    i,j = getBoundsIndices(z, z0)

    # Interpolate for f(x0,y0) holding z fixed
    fi = interpolate2D(f[:,:,z[i]], x, y, x0, y0)
    fj = interpolate2D(f[:,:,z[j]], x, y, x0, y0)

    # Interpolate between the planes
    return interpolate1D([fi, fj], [z[i], z[j]], z0)

end

function getBoundsIndices(x::Vector, x0)::Tuple
    """Returns indices i and j of points such that x[i] < x <= x[j].
    Values of x0 > max are treated periodically. Assumes x is ascending."""

    # Make sure x0 is in range
    if minimum(x) <= x0 < maximum(x)
        @debug "Point out of range of input vectors"
        x0 = mod1(x0, maximum(x))
    end

    # Get upper bounds
    j = findfirst(x .>= x0)
    i = j == 1 ? length(x) : j-1

    return (i, j)
end


using Plots
function interpolate2DTest()
    
    data = rand(Float64, (4,4))
    x = Vector(1:4)
    y = Vector(1:4)

    xpoints = 1:0.1:4
    ypoints = 1:0.1:4

    int2d = [interpolate2D(data, x, y, x0, y0) for x0 in xpoints, y0 in ypoints]
    display(surface(xpoints, ypoints, int2d))
end

function interpolate3DTest()
    
    data = rand(Float64, (4,4,4))
    x = Vector(1:4)
    y = Vector(1:4)
    z = Vector(1:4)

    xpoints = 1:0.1:4
    ypoints = 1:0.1:4
    zpoints = 1:0.1:4

    [interpolate3D(data, x, y, z, x0, y0, z0) for x0 in xpoints, y0 in ypoints, z0 in zpoints]
end

interpolate2DTest()
interpolate3DTest()