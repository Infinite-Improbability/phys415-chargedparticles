using LinearAlgebra: cross, dot, norm
using Plots
# add unitful?

mutable struct Particle
    r::Vector{Real} # position
    v::Vector{Real} # velocity
    q::Real
    m::Real
end

# Define electric field as function of position
E(r::Vector{Real})::Vector{Real} = [0, 0, 0.01]

# Define magnetic field as function of position
B(r::Vector{Real})::Vector{Real} = [0, 0, 0.01]

function stepVelocity!(part::Particle, dt::Real)
    halfEV = part.q / part.m * E(part.r) * dt / 2 # half electric acceleration

    t = part.q / part.m * B(part.r) * dt / 2
    s = 2 * t / (1 + dot(t,t)) # s= 2t/(1+|t|^2)

    vMinus = part.v + halfEV
    vPrime = vMinus + cross(vMinus, t)
    vPlus = vMinus + cross(vPrime, s)
    part.v = vPlus + halfEV
end

function stepPosition!(part::Particle, dt::Real)
    part.r += part.v * dt
end

function lamorRadius(part::Particle)::Real
    # First want the tangential velocity
    n = part.r / norm(part.r) # unit vector of position, normal to surface
    vTang = norm(part.v - dot(part.v, n) * n) # v⟂ = v - v∥ = v - (v⋅n)n
    
    return part.m * vTang / part.q / norm(B(part.r))
end

# Init system
dt = 3e-11
iterations = 1000
p1 = Particle([0,0,0], [0,1e5,0], -1.602e-19, 9.109e-31)
stepVelocity!(p1, -dt/2) # Offset velocity back 1/2 step

# Init plot
# We'll base the limits on the Larmor radius of the particle
lim = 1.1 * lamorRadius(p1)
plt = plot3d(
    1, # a single empty series
    xlim = (-lim, lim),
    ylim = (-lim, lim),
    zlim = (-lim, lim),
    title = "Charged particle goes zoom",
    marker = 2,
    legend=:none
)

# Main loop
anim = @gif for i in 1:iterations
    stepVelocity!(p1, dt)
    stepPosition!(p1, dt)
    push!(plt, p1.r...) # expands out positions
end every 10