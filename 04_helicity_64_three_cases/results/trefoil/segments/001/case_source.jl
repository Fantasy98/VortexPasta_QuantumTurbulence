#!/usr/bin/env julia
# Installed TrefoilKnot has coordinates (sin(s)+2sin(2s), cos(s)-2cos(2s), -sin(3s)).
# Run: julia --project=. --threads=2 trefoil.jl
include("simulation_common.jl")

function initial_curves()
    knot = define_curve(TrefoilKnot(); scale=RADIUS/3,
        translate=(π, π, π), orientation=1)
    [arclength_filament(knot, 128)]
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_case("trefoil", initial_curves, @__FILE__)
end
