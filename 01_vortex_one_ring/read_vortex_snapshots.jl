# Used by render_vortex_ring.py to read existing VTKHDF files without GLMakie.
# Each output line is a JSON array of filaments, each containing [x, y, z] points.
# Uses only the VortexPasta dependency already in this project.
using VortexPasta.FilamentIO: read_vtkhdf
using VortexPasta.Filaments: QuinticSplineMethod

for filename in ARGS
    filaments = read_vtkhdf(filename, Float64, QuinticSplineMethod())
    print('[')
    for (filament_id, filament) in enumerate(filaments)
        filament_id > 1 && print(',')
        print('[')
        # Include the closing endpoint, as in the original VTKHDF geometry.
        for i in 1:(length(filament) + 1)
            i > 1 && print(',')
            point = filament[i]
            all(isfinite, point) || error("Non-finite coordinates in $filename")
            print('[', join(point, ','), ']')
        end
        print(']')
    end
    println(']')
end
