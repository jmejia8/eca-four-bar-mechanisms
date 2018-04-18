include("findStructure.jl")

function example()
    ncase = 3

    println("Optimize...")
    p, er, precisionpts = findStructure(ncase)


    println("Generate animation...")
    img = animateFourBarMechanism(p;
                    precision_points=precisionpts,
                    title   = @sprintf("Error: %e", er),
                    xlimits = (-5, 5),
                    ylimits = (-5, 5)
                    )

    run(`xviewer $(img.filename)`)

end

example()

