using Javis
using Unitful

function makeRender(iterations::Int, particles::Vector, positions::Array{Quantity,3}; bgcolour="white", lenUnit::Unitful.LengthFreeUnits=u"m", outputFile="circle.gif", size::Int=500)
    print("Beginning animated render.\n")
    # Init video
    function ground(args...)
        background(bgcolour) # canvas background
        sethue("black") # pen color
    end

    function object(p=O, colour="black")
        sethue(colour)
        circle(p, 5, :fill)
        return p
    end

    scaleFactor = size / ustrip(lenUnit, maximum(positions[:,:,:2]))

    myvideo = Video(size, size)
    Background(1:iterations, ground)

    print("Render initialisation complete. Loading points.\n")

    # Get paths followed by particles [iterations, particles]
    point_positions = Matrix{Point}(undef, iterations, length(particles))
    for j in 1:length(particles)
        for i in 1:iterations
            p = ustrip.(lenUnit, positions[i, j, :]) * scaleFactor
            point_positions[i, j] = Point(p[1], p[2])
        end
    end

    # Add particles to vide
    dots = Vector{Object}(undef, length(particles))
    for i in 1:length(particles)
        p = particles[i]
        colour = deep[ustrip(u"C", p.q)]
        pos = ustrip.(lenUnit, p.r) * scaleFactor
        dots[i] = Object((args...) -> object(O, colour), Point(pos[1], pos[2]))
        act!(dots[i], Action(1:iterations, follow_path(point_positions[:, i])))
    end

    print("Points added to video. Beginning render.\n")

    # Render output
    render(
        myvideo;
        pathname=outputFile
    )

    print("Render complete. File saved as $(outputFile).\n")
end