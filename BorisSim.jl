using LinearAlgebra: cross, dot, norm
using Unitful
using Random
using GLMakie
import ColorSchemes.deep
include("interpolation.jl")

mutable struct Particle
    """Properties of a single particle in the system."""
    r::Vector # position vector [x, y, z]
    v::Vector # velocity vector
    q::Quantity # charge
    m::Quantity # mass
end

# Allow copying particle structs
Base.copy(s::Particle) = Particle(s.r, s.v, s.q, s.m)
# Random particle generation
Random.rand(rng::AbstractRNG, ::Random.SamplerType{Particle}) = Particle([0u"m", 0u"m", 0u"m"],
rand(rng, typeof(1.0u"m/s"), 3), rand(rng, typeof(1.0u"C")), rand(rng, typeof(1.0u"kg")))

# Define the electric and magnetic fields.
# They are always defined as functions of position - even for
# constant fields - so that you don't have to go through changing
# the Boris method code whenever you switch to a position dependent field.
# The field functions can be anything that returns a three element vector with units
# appropriate to the field type.
# Units are denoted using Unitful's notation.
# e.g. 1u"V/m" is 1 Volt per metre
E(r::Vector)::Vector = [0u"V/m", 0u"V/m", 0.01u"V/m"]
B(r::Vector)::Vector = [0u"T", 0u"T", 1000000000u"T"]

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
    part.r = mod1.(ustrip.(u"m", part.r), lim) * 1u"m"
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
iterations = 10000

# Define particles, in inital state
particles = rand(Particle, 4)
# Set lims where everything loops around
lim = 10

# Init plot
GLMakie.activate!()
scene = Scene(;
    # clear everything behind scene
    clear = true,
    # the camera struct of the scene.
    visible = true,
    # ssao and light are explained in more detail in `Documetation/Lighting`
    ssao = Makie.SSAO(),
    # Creates lights from theme, which right now defaults to `
    # set_theme!(lightposition=:eyeposition, ambient=RGBf(0.5, 0.5, 0.5))`
    lights = Makie.automatic,
    backgroundcolor = :white,
    resolution = (500, 500)
)
cam3d!(scene)

spheres = [Sphere(Point3f(ustrip.(u"m", p.r)...), ustrip(u"kg", p.m)) for p in particles]
colours = deep[[ustrip(u"C", p.q) for p in particles]]
sphere_plots = [mesh!(scene, spheres[i], color=colours[i]) for i in 1:length(particles)]

# Offset initial velocity back 1/2 step so it leapfrogs with position
stepVelocity!.(particles, -dt / 2)

# Main loop
print("Beginning loop.\n")
record(scene, "scene.mp4", 1:iterations) do i
    # ! indicates inplace functions
    stepVelocity!.(particles, dt)
    stepPosition!.(particles, dt)
    px = [ustrip.(u"m", p.r[1]) for p in particles]
    py = [ustrip.(u"m", p.r[2]) for p in particles]
    pz = [ustrip.(u"m", p.r[3]) for p in particles]
    translation_vectors = Vec3f.(px, py, pz)
    translate!.(sphere_plots, translation_vectors)
end

print("Loop complete.\n")
