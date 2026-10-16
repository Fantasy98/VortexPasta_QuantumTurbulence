#!/usr/bin/env julia
# Two perpendicular PredefinedCurves.Ring circles; oriented linking number -1.
# Run: julia --project=. --threads=2 hopf_link.jl
include("simulation_common.jl")

function initial_curves()
    first = define_curve(Ring(); scale=RADIUS,
        translate=(π-RADIUS/2, π, π), orientation=1)
    # Rotate +π/2 around x: (cos(s), sin(s), 0) -> (cos(s), 0, sin(s)).
    rotate_x = [1.0 0.0 0.0; 0.0 0.0 -1.0; 0.0 1.0 0.0]
    second = define_curve(Ring(); scale=RADIUS, rotate=rotate_x,
        translate=(π+RADIUS/2, π, π), orientation=1)
    [arclength_filament(first, 64), arclength_filament(second, 64)]
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_case("hopf_link", initial_curves, @__FILE__)
end
