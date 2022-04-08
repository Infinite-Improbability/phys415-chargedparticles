using LinearAlgebra: cross, dot, norm
using Plots
using Unitful
using UnitfulRecipes

mutable struct Particle
    """Properties of a single particle in the system."""
    r::Vector # position vector [x, y, z]
    v::Vector # velocity vector
    q::Quantity # charge
    m::Quantity # mass
end

# Define electric field as function of position
# E(r::Vector)::Vector = [0u"V/m", 0u"V/m", 0.01u"V/m"]
cutoff = -4e-7u"m"
E(r::Vector)::Vector = r[3] > cutoff ? [0u"V/m", 0u"V/m", 0.01u"V/m"] : [0u"V/m", 0u"V/m", 0.00u"V/m"]
# We aren't actually using the position here but by passing it we make
# it easier in the future

# Define magnetic field as function of position
# B(r::Vector)::Vector = [0u"T", 0u"T", 0.01u"T"]
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
    n = B(p1.r) / norm(B(p1.r)) # unit vector of magnetic field
    vTang = norm(part.v - dot(part.v, n) * n) # v⟂ = v - v∥ = v - (v⋅n)n
    return part.m * vTang / part.q / norm(B(part.r))
end

# Init system
dt = 3e-11u"s" # timestep
iterations = 1000
# p1 = Particle([0u"m", 0u"m", 0u"m"], [0u"m/s", 1e5u"m/s", 0u"m/s"], -1.602e-19u"C", 9.109e-31u"kg")
p1 = Particle([1.0u"m", 0u"m", 0u"m"], [0.0u"m/s", 1u"m/s", 0u"m/s"], -1.602e-19u"C", 9.109e-31u"kg")
p1.r[1] = -lamorRadius(p1) # adjust particle position so that it rotates around origin. Polarity dependent
display(p1)

# Analytical solutions
time = 0u"s" : dt : (dt*(iterations-1))

# E=0, B=const, assumes orbit around origin, B up in z only
ω = p1.q * norm(B(p1.r)) / p1.m # gyrofrequency
pos = zeros(typeof(1.0u"m"), (1000,3))
for i in 1:iterations
    θ = -ω*time[i]
    pos[i,:] = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] * p1.r # rotate about z
end

# E=[0,0,0.01] when z>cutoff, B=0
accel = p1.q * E([0.0u"m", 0u"m", 0u"m"]) / p1.m
tField = sqrt(abs(2*cutoff/accel[3])) # time to clear E field
vField = accel*tField # velocity gained in this time
Δd(t) = t <= tField ? (t*p1.v+0.5*accel*t^2)[3] : cutoff + vField[3]*(t-tField)
pos1 = (time*p1.v')[:,3] + Δd.(time) # row i is [x y z] at time[i]

# Combine E and B contributions
pos[:,3] += pos1

# Init plot
# We'll base the limits on the Larmor radius of the particle
# lim = ustrip(1.1 * lamorRadius(p1)) # Fails if initial pos is origin, code executes with auto limits
plt = plot3d(
    1, # a single empty series
    # lims=(-lim, lim),
    # title="E=const B=0",
    marker=2,
    legend=:none,
    xlabel="x",
    ylabel="y",
    zlabel="z",
    #line_z=ustrip.(time), # colour line based on time
    #colorbar=:bottom,
    label="Numerical",
    #markershape=:none
)

# Offset initial velocity back 1/2 step so it leapfrogs with position
stepVelocity!(p1, -dt / 2) 

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

# The gif macro doesn't allow a custom output location.
# It would be possible if I went more hands on in
# the gif construction but it doesn't seem worth the effort.
# You cqn export the last frame as an image with
# savefig("filename")
print("GIF saved at $(anim.filename)\n")
display(anim)

#Add analytical solution to plot
plot!(pos[:,1], pos[:,2], pos[:,3],
    label="Analytical",
    legend=:right,
    # line_z=ustrip.(time),
    lims=:auto)
display(plt)

print("Magnitude of error in final positions: $(norm(pos[end,:]-p1.r))\n")