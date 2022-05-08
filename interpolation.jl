function interpolate1D(f::Vector, x::Vector, x0)
    """Returns value of function f at x0 by linear interpolation of provided datapoints.
    x should be ascending. Interpolation is periodic"""

    # Get an upper bound index for x0
    # We use mod1 to make it periodic
    i = mod1(findfirst(x .>= x0), length(x))

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

# using Plots
# function interpolate2DTest()
    
#     data = rand(Float64, (4,4))
#     x = Vector(1:4)
#     y = Vector(1:4)

#     xpoints = 1:0.1:4
#     ypoints = 1:0.1:4

#     int2d = [interpolate2D(data, x, y, x0, y0) for x0 in xpoints, y0 in ypoints]
#     surface(xpoints, ypoints, int2d)
# end

# interpolate2DTest()