using LinearAlgebra: cross, dot, norm
using Plots
using Unitful

mutable struct Particle
    """Properties of a single particle in the system."""
    r::Vector # position vector [x, y, z]
    v::Vector # velocity vector
    q::Quantity # charge
    m::Quantity # mass
end

# Define electric field as function of position
E(r::Vector)::Vector = [0u"V/m", 0u"V/m", 0.01u"V/m"]
# We aren't actually using the position here but by passing it we make
# it easier in the future

# Define magnetic field as function of position
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
    n = part.r / norm(part.r) # unit vector of position, normal to surface
    vTang = norm(part.v - dot(part.v, n) * n) # v⟂ = v - v∥ = v - (v⋅n)n

    return part.m * vTang / part.q / norm(B(part.r))
end

# Init system
dt = 3e-11u"s" # timestep
iterations = 1000
p1 = Particle([0u"m", 0u"m", 0u"m"], [0u"m/s", 1e5u"m/s", 0u"m/s"], -1.602e-19u"C", 9.109e-31u"kg")
stepVelocity!(p1, -dt / 2) # Offset initial velocity back 1/2 step

# Init plot
# We'll base the limits on the Larmor radius of the particle
lim = 1.1 * lamorRadius(p1)
plt = plot3d(
    1, # a single empty series
    xlim=(-lim, lim),
    ylim=(-lim, lim),
    zlim=(-lim, lim),
    title="Charged particle goes zoom",
    marker=2,
    legend=:none,
    xlabel="x (m)",
    ylabel="y (m)",
    zlabel="z (m)"
)

# Main loop
anim = @gif for i in 1:iterations
    # ! indicates inplace functions
    stepVelocity!(p1, dt)
    stepPosition!(p1, dt)
    push!(plt, ustrip.(u"m", p1.r)...) # expands out positions
    # Pushing Quantities to the plot failed even when using UnitfulRecipes
    # (a package to give Plots unit support)
    # This at least works.
end every 10 # save every tenth run as a frame

print("GIF saved at $(anim.filename)")
display(anim)
# The gif macro doesn't allow a custom output location.
# It would be possible if I went more hands on in
# the gif construction but it doesn't seem worth the effort.