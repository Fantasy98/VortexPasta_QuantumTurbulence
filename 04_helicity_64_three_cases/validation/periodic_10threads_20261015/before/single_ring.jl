#!/usr/bin/env julia
# Run: julia --project=. --threads=2 single_ring.jl
include("simulation_common.jl")

function initial_curves()
    ring = define_curve(Ring(); scale=RADIUS, translate=(π, π, π), orientation=1)
    [arclength_filament(ring, 64)]
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_case("single_ring", initial_curves, @__FILE__)
end
