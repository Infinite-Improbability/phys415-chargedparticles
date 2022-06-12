using LinearAlgebra: cross, dot, norm
using Unitful
using Random
include("readData.jl") # Read field in from file
include("renders.jl") # Functions to generate outputs
include("interpolation.jl")

# -----------------------------------------------------------------------
mutable struct Particle
    """Properties of a single particle in the system."""
    r::Vector # position vector [x, y, z]
    v::Vector # velocity vector
    q # charge
    m # mass
end

# Create a way to generate random particles
Random.rand(rng::AbstractRNG, ::Random.SamplerType{Particle}) = Particle(
    (rand(rng, typeof(1.0u"m"), 3) - rand(rng, typeof(1.0u"m"), 3)) * 1e-7,
    rand(rng, typeof(1.0u"m/s"), 3) - rand(rng, typeof(1.0u"m/s"), 3),
    rand(rng, typeof(1.0u"C")) - rand(rng, typeof(1.0u"C")),
    rand(rng, typeof(1.0u"kg"))
)

# More limited random values, makes cooler plots
# Random.rand(rng::AbstractRNG, ::Random.SamplerType{Particle}) = Particle(
#     [0u"m", 0u"m", 0u"m"],
#     rand(rng, typeof(1.0u"m/s"), 3),
#     rand(rng, typeof(1.0u"C")),
#     rand(rng, typeof(1.0u"kg"))
# )


# -----------------------------------------------------------------------

# Define the electric and magnetic fields.
# They are always defined as functions of position - even for
# constant fields - so that you don't have to go through changing
# the Boris method code whenever you switch to a position dependent field.

# The field functions can be anything of the form E(r::Vector)::Vector{T}
# The output is a vector of length 3 and in units apprioate to the field type
# There must be an E function for the electric field and a B function for the magnetic field

# Units are denoted using Unitful's notation.
# e.g. 1u"V/m" is 1 Volt per metre
# Just uncomment whichever block you want to use.

# Uniform fields
E(r::Vector)::Vector = [0u"V/m", 0u"V/m", 0.01u"V/m"]
B(r::Vector)::Vector = [0u"T", 0u"T", 1000000000u"T"]

# Load magnetic field in from HDF5 data
# Interpolation is slow so don't use it with lots of particles
# E(r::Vector)::Vector = [0u"V/m", 0u"V/m", 0.0u"V/m"]
# Bx, By, Bz = loadHDF5()
# # We'll make up the spatial dimensions, in metres
# x = range(-128, 128, 256) * 1e-8
# y = range(-128, 128, 256) * 1e-8
# z = range(-128, 128, 256) * 1e-8
# function B(r::Vector)::Vector
#     r = ustrip.(u"m", r)
#     fx = interpolate3D(Bx, x, y, z, r...)
#     fy = interpolate3D(By, x, y, z, r...)
#     fz = interpolate3D(Bz, x, y, z, r...)
#     return [fx, fy, fz] * 100000u"T" # Scale factor was chosen to get clear behaviour for default particle generation
# end

# Interpolate over a uniform field
# E(r::Vector)::Vector = [0u"V/m", 0u"V/m", 0.0u"V/m"]
# Bx = fill(100, (256,256,256))
# By = fill(100, (256,256,256))
# Bz = fill(100, (256,256,256))
# # We'll make up the spatial dimensions, in metres
# x = range(-128, 128, 256) * 1e-7
# y = range(-128, 128, 256) * 1e-7
# z = range(-128, 128, 256) * 1e-7
# function B(r::Vector)::Vector
#     r = ustrip.(u"m", r)
#     fx = interpolate3D(Bx, x, y, z, r...)
#     fy = interpolate3D(By, x, y, z, r...)
#     fz = interpolate3D(Bz, x, y, z, r...)
#     return [fx, fy, fz] * 100000u"T" # Scale factor was chosen to get clear behaviour for default particle generation
# end


# -----------------------------------------------------------------------

function stepVelocity!(part::Particle, dt::Quantity)
    """Update velocity with change over timestep dt using Boris method."""
    halfEV = part.q / part.m * E(part.r) * dt / 2 # half electric acceleration

    t = part.q / part.m * B(part.r) * dt / 2
    s = 2 * t / (1 + dot(t, t)) # s= 2t/(1+|t|^2)

    vMinus = part.v + halfEV
    vPrime = vMinus + cross(vMinus, t)
    vPlus = vMinus + cross(vPrime, s)
    part.v = vPlus + halfEV
end


function stepPosition!(part::Particle, dt::Quantity)
    """Update position over timestep dt with current velocity."""
    part.r += part.v * dt
end


function lamorRadius(part::Particle)::Quantity
    """Calculate lamorRadius of particle."""
    # First want the tangential velocity
    n = B(part.r) / norm(B(part.r)) # unit vector of magnetic field
    vTang = norm(part.v - dot(part.v, n) * n) # v⟂ = v - v∥ = v - (v⋅n)n
    return part.m * vTang / part.q / norm(B(part.r))
end

# -----------------------------------------------------------------------

# Init system
print("Initialising system.\n")
dt = 3e-11u"s" # timestep
iterations = 10000
# Define particles, in inital state
# The second argument to rand sets the number of particles.
particles = rand(Particle, 100)
print("There are $(length(particles)) particles.\n")
# Offset initial velocity back 1/2 step so it leapfrogs with position
stepVelocity!.(particles, -dt / 2)
# We'll store history in the following 3D array with indices [time, particle, coordinate]
positions = Array{Quantity,3}(undef, iterations, length(particles), 3)

# Main loop
print("Init complete. Beginning loop.\n")
for i in 1:iterations
    stepVelocity!.(particles, dt)
    stepPosition!.(particles, dt)
    positions[i, :, :] = [p.r[j] for p in particles, j in 1:3]
end
print("Loop complete.\n")

print("Making graphics.\n")
plotTrajectories(positions)
print("Close plot and press enter to load next plot.")
readline()
plotAtTime(positions, iterations, dt)
print("Close plot and press enter to load next plot.")
readline()
plotParticle(positions, 1)

# This is buggy but the output is the most interesting of the various options
# makeRender(iterations, particles, positions, size=100)

print("Graphics complete. Press enter to exit.")
readline() # Doesn't interact well with running in REPL.