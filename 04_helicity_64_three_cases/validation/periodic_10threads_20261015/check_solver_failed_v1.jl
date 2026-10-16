# Short validation of the existing three configurations, not new production runs.
include(joinpath(@__DIR__, "..", "..", "simulation_common.jl"))

function solver_for(fs, p, cfg)
    hmin, hmax = cfg["hmin"], cfg["hmax"]
    adapt = AdaptBasedOnVelocity(0.15min(DREC,hmin)) | MaximumTimestep(DT_MAX) |
        AdaptBasedOnSegmentLength(0.25)
    init(VortexFilamentProblem(fs, 0.1, p), RK4(); dt=DT_MAX,
        refinement=RefineBasedOnSegmentLength(hmin,hmax),
        reconnect=ReconnectFast(DREC;nthreads=1,max_passes=2,use_velocity=true),
        adaptivity=adapt, fold_periodic=false)
end

function main()
    @assert Threads.nthreads() == 10
    @assert BLAS.get_num_threads() == 1
    p = parameters()
    @assert p.Ls == (BOX,BOX,BOX)
    @assert p.longrange.Ns == (64,64,64)
    report = Dict{String,Any}("Julia_threads"=>Threads.nthreads(),
        "BLAS_threads"=>BLAS.get_num_threads(), "periodic_Ls"=>collect(p.Ls),
        "mesh"=>collect(p.longrange.Ns), "julia"=>string(VERSION),
        "VortexPasta"=>string(pkgversion(VortexPasta)))
    started = time()
    for name in ("single_ring", "hopf_link", "trefoil")
        dir = joinpath(ROOT,"results",name)
        cfg = TOML.parsefile(joinpath(dir,"config.toml"))
        checkpoint = load_checkpoint(joinpath(dir,"segments","001","checkpoint_00000000.vtkhdf"),
            Float64,QuinticSplineMethod())
        iter = solver_for(deepcopy(checkpoint.fs),p,cfg)
        row = Dict{String,Any}("initial_nodes"=>length.(iter.fs))
        # A whole-cell translation through either z face must preserve velocity.
        errors = Float64[]
        for sign in (-1,1)
            shifted = deepcopy(checkpoint.fs)
            for f in shifted
                for i in eachindex(f)
                    f[i] += Vec3(0.,0.,sign*BOX)
                end
                Filaments.update_coefficients!(f)
            end
            image_iter = solver_for(shifted,p,cfg)
            error = maximum(norm(iter.vs[j][i]-image_iter.vs[j][i])
                for j in eachindex(iter.fs) for i in eachindex(iter.fs[j]))
            scale = maximum(norm(v) for vs in iter.vs for v in vs)
            push!(errors,error/scale)
            @assert error/scale < 1e-9
        end
        row["z_minus_plus_period_velocity_relative_errors"] = errors
        integration_start=time()
        for _ in 1:20
            step!(iter)
        end
        row["integration_elapsed_seconds"] = time()-integration_start
        @assert iter.nstep == 20 && iter.t > 0
        @assert all(isfinite, diagnostic(iter))
        row["steps"] = iter.nstep
        row["final_time"] = iter.t
        row["final_nodes"] = length.(iter.fs)
        row["finite_diagnostics"] = true
        report[name] = row
        println(name, ": 20 steps on ", Threads.nthreads(), " threads, t=", iter.t,
            "; periodic image relative velocity errors=", errors)
        flush(stdout)
    end
    report["validation_elapsed_seconds"] = time()-started
    report["note"] = "Short ten-thread checks only; production t=10 results are unchanged. No speedup benchmark."
    file=joinpath(@__DIR__,"solver_checks.toml")
    isfile(file) && error("Preserving existing validation results")
    save_toml(file,report)
end
main()
