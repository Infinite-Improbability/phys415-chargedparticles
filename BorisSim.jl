using LinearAlgebra: cross, dot, norm
using Plots
using Unitful
using UnitfulRecipes
include("interpolation.jl")

mutable struct Particle
    """Properties of a single particle in the system."""
    r::Vector # position vector [x, y, z]
    v::Vector # velocity vector
    q::Quantity # charge
    m::Quantity # mass
end

# Define the electric and magnetic fields.
# They are always defined as functions of position - even for
# constant fields - so that you don't have to go through changing
# the Boris method code whenever you switch to a position dependent field.
# The field functions can be anything that returns a three element vector with units
# appropriate to the field type.
# Units are denoted using Unitful's notation.
# e.g. 1u"V/m" is 1 Volt per metre
E(r::Vector)::Vector = [0u"V/m", 0u"V/m", 0.01u"V/m"]
B(r::Vector)::Vector = [0u"T", 0u"T", 0.01u"T"]

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


# Init system
print("Initialising system and plot.\n")
dt = 3e-11u"s" # timestep
iterations = 1000

# Define particle, in inital state.
p1 = Particle(
    [0u"m", 0u"m", 0u"m"], # Initial position
    [0u"m/s", 1e5u"m/s", 0u"m/s"], # Inital velocity
    -1.602e-19u"C", # Charge
    9.109e-31u"kg" # Mass
)

# Init plot
# Plot attributes have been chosen to cleanly represent
# a single particle with no analytical line.
# For more options see https://docs.juliaplots.org/stable/attributes/
plt = plot3d(
    1, # a single empty series
    # title="Some title",
    legend=:none,
    xlabel="x (m)",
    ylabel="y (m)",
    zlabel="z (m)",
    label="Numerical",
    markershape=:none
)

# Offset initial velocity back 1/2 step so it leapfrogs with position
stepVelocity!(p1, -dt / 2)

# Main loop
print("Beginning loop.\n")
anim = @gif for i in 1:iterations
    # ! indicates inplace functions
    stepVelocity!(p1, dt)
    stepPosition!(p1, dt)
    push!(plt, ustrip.(u"m", p1.r)...) # ... expands out positions. ustrip applied elementwise due to . before brackets
    # Pushing Quantities to the plot failed even when using UnitfulRecipes
    # (a package to give Plots unit support)
    # This at least works.
end every 10 # save every tenth run as a frame
print("Loop complete.\n")

# The gif macro doesn't allow a custom output location.
# It would be possible if I went more hands on in
# the gif construction but it doesn't seem worth the effort.
# You cqn export the last frame as an image with
# savefig("filename")
print("GIF saved at $(anim.filename)\n")
# This is probably unneeded, as Julia Plots prints this anyway,
# But I include it just in case.

# The GUI display can only handle the last frame, not the animation.
# It also doesn't seem to work when running in the REPL from VS Code,
# but is helpful if funning in the terminal.
gui()
print("Press any key to exit ") # Don't immediately close on us
readline()