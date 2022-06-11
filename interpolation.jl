function interpolate1D(f::Vector, x::Vector, x0)
    """Returns value of function f at x0 by linear interpolation of provided datapoints.
    x should be ascending. Interpolation is periodic"""

    # Get an upper bound index for x0
    # We use mod1 to make it periodic
    i = max(mod1(findfirst(x .>= x0), length(x)), 2)

    # Interpolate and return
    f[i-1] + (f[i] - f[i-1])*(x0 - x[i-1])/(x[i] - x[i-1])
end

function interpolate2D(f::Matrix, x::Vector, y::Vector, x0, y0)
    """Returns value of function f at (x0,y0) by bilinear interpolation of provided datapoints.
    x should be ascending. Interpolation is periodic."""

    # Find lower indices of grid square containing (x0, y0)
    i = mod1(findfirst(x .>= x0), length(x) - 1)
    j = mod1(findfirst(y .>= y0), length(y) - 1)

    # Prepare upper bounds so we don't go out of bounds
    ip = mod1(i+1, length(x))
    jp = mod1(j+1, length(y))

    t = (x0 - x[i])/(x[ip] - x[i])
    u = (y0 - y[i])/(y[ip] - y[i])

    t1 = 1-t
    u1 = 1-u

    # Interpolate and return
    t1*u1*f[i,j] + t*u1*f[ip,j] + t*u*f[ip,jp] + t1*u*f[i,jp]
end

function interpolate2DV2(f::Matrix, x::Vector, y::Vector, x0, y0)
    """
    Returns value of function f at (x0,y0) by bilinear interpolation of provided datapoints.
    x should ascend as row index increases. y should ascend as column index increases.
    Interpolation is periodic.
    """

    # Find upper and lower y indices of grid square containing (x0, y0)
    jp = mod1(findfirst(y .>= y0), length(y))
    j = jp-1

    # Intepolate for x at borders of grid (fixed y)
    fj = interpolate1D(f[:,j], x, x0)
    fjp = interpolate1D(f[:,jp], x, x0)

    return interpolate1D([fj, fjp], [y[j], y[jp]], y0)
end

function interpolate3D(f::Matrix, x::Vector, x0, y0, z0)
    """Returns value of function f at (x0,y0,z0) by bilinear interpolation of provided datapoints.
    x should be ascending. Interpolation is periodic."""

    # Find lower indices of grid square containing (x0, y0)
    i = mod1(findfirst(x .>= x0), length(x) - 1)
    j = mod1(findfirst(y .>= y0), length(y) - 1)
    k = mod1(findfirst(z .>= z0), length(z) - 1)

    # Prepare upper bounds so we don't go out of bounds
    ip = mod1(i+1, length(x))
    jp = mod1(j+1, length(y))
    kp = mod1(k+1, length(k))

    t = (x0 - x[i])/(x[ip] - x[i])
    u = (y0 - y[i])/(y[ip] - y[i])
    v = (z0 - z[i])/(z[ip] - z[ip])

    t1 = 1-t
    u1 = 1-u
    z1 = 1-v
end

function getBoundsIndices(x::Vector, x0)::Tuple
    """Returns indices i and j of points such that x[i] < x <= x[j].
    Values of x0 > max are treated periodically. Assumes x is ascending."""

    # Make sure x0 is in range
    if min(x) <= x0 < max(x)
        @debug "Point out of range of input vectors"
    end
    x0 = mod1(x0, max(x))

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

    int2dv2 = [interpolate2DV2(data, x, y, x0, y0) for x0 in xpoints, y0 in ypoints]
    display(surface(xpoints, ypoints, int2dv2))
end

interpolate2DTest()