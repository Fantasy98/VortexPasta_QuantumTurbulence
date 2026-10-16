# Shared numerical settings and output for the three analytic base cases.
using VortexPasta
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.PredefinedCurves: Ring, TrefoilKnot, define_curve
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using LinearAlgebra, Statistics, Printf, TOML, SHA, Dates

BLAS.set_num_threads(1)
const ROOT = @__DIR__
const COMMON_SOURCE = read(@__FILE__, String)
const RADIUS = 1.0
const BOX = 2π
const GAMMA = 4π*2*0.1171875/sqrt(2)
const CORE = 0.001
const DELTA = 0.615
const END_TIME = 10.0
const MESH = (64, 64, 64)  # Fourier/Ewald mesh for ALL three cases, including trefoil.
const BETA = 4.5
const QUAD_ORDER = 10
const DREC = 0.08
const DT_MAX = 0.001
const DIAG_INTERVAL = 0.02
const SNAPSHOT_INTERVAL = 0.05
const CHECKPOINT_INTERVAL = 1.0

"""Sample a PredefinedCurves function approximately uniformly in arclength.

Only its parameter is remapped; the positions still come from the actual
PredefinedCurves function, with no imported centerlines or replacement curve.
"""
function arclength_filament(curve, n)
    ts = collect(range(0.0, 1.0; length=16385))
    xs = curve.(ts)
    distances = [0.0; cumsum(norm(xs[i+1]-xs[i]) for i in 1:length(xs)-1)]
    function arclength_curve(t)
        d = clamp(t, 0.0, 1.0)*last(distances)
        j = clamp(searchsortedlast(distances, d), 1, length(ts)-1)
        w = (d-distances[j])/(distances[j+1]-distances[j])
        curve((1-w)*ts[j]+w*ts[j+1])
    end
    Filaments.init(arclength_curve, ClosedFilament{Float64}, n, QuinticSplineMethod())
end

function parameters()
    alpha = (MESH[1]-1)/(4BETA)
    # Ls sets periodic physics on all axes, including matching z=0 and z=2π.
    ParamsBiotSavart(;
        Γ=GAMMA, a=CORE, Δ=DELTA, Ls=(BOX,BOX,BOX), Ns=MESH,
        α=alpha, rcut=BETA/alpha,
        quadrature=GaussLegendre(QUAD_ORDER),
        quadrature_near_singularity=GaussLegendre(QUAD_ORDER),
        lia_segment_fraction=0.2,
        backend_short=CellListsBackend(2),
        backend_long=NonuniformFFTsBackend(CPU(); σ=1.5, m=HalfSupport(6)))
end

function write_curves(path, fs)
    open(path, "w") do io
        println(io, "loop,node,x,y,z")
        for (j,f) in enumerate(fs), i in eachindex(f)
            println(io, join((j,i,f[i]...), ','))
        end
    end
end

function diagnostic(iter)
    fs=iter.fs
    ds=[norm(f[i+1]-f[i]) for f in fs for i in eachindex(f)]
    curvature=[f[i,CurvatureScalar()] for f in fs for i in eachindex(f)]
    c=isempty(fs) ? Vec3(NaN,NaN,NaN) : mean(x for f in fs for x in nodes(f))
    impulse=isempty(fs) ? Vec3(0.,0.,0.) : Diagnostics.vortex_impulse(iter)
    [iter.t, GAMMA*iter.t/RADIUS^2, iter.nstep, iter.dt, length(fs), sum(length,fs;init=0),
     isempty(ds) ? NaN : minimum(ds), isempty(ds) ? NaN : mean(ds), isempty(ds) ? NaN : maximum(ds),
     isempty(curvature) ? NaN : maximum(curvature),
     isempty(fs) ? 0. : Diagnostics.filament_length(iter),
     isempty(fs) ? 0. : Diagnostics.kinetic_energy(iter), impulse...,
     isempty(fs) ? 0. : Diagnostics.helicity(iter)/GAMMA^2, c...,
     iter.stats.reconnection_count, iter.stats.reconnection_length_loss,
     iter.stats.refinement_length_loss, iter.stats.filaments_removed_count,
     iter.stats.filaments_removed_length]
end

function save_toml(path, value)
    open(path,"w") do io; TOML.print(io,value); end
end

function run_case(name, make_curves, entry_file)
    all(x in ("--check", "--resume") for x in ARGS) || error("Options: --check or --resume")
    fs0=make_curves()
    ds0=[norm(f[i+1]-f[i]) for f in fs0 for i in eachindex(f)]
    hmin=0.65minimum(ds0); hmax=1.4minimum(ds0)
    if "--check" in ARGS
        @assert all(isfinite, (x[k] for f in fs0 for x in nodes(f) for k=1:3))
        println(name, ": ", length(fs0), " loops, nodes=", length.(fs0),
            ", segments=", extrema(ds0), ", mesh=", MESH, ", t_end=", END_TIME,
            ", periodic_Ls=", (BOX,BOX,BOX), ", Julia_threads=", Threads.nthreads(),
            ", BLAS_threads=", BLAS.get_num_threads(), ", VortexPasta=", pkgversion(VortexPasta))
        return
    end
    dir=joinpath(ROOT,"results",name)
    resume="--resume" in ARGS
    if !resume && isdir(dir) && !isempty(readdir(dir))
        error("Existing results preserved: $dir. Use --resume only for an unfinished case.")
    end
    mkpath(joinpath(dir,"segments"))
    segments=sort(filter(s->isdir(joinpath(dir,"segments",s)),readdir(joinpath(dir,"segments"))))
    restart=nothing
    if resume
        isempty(segments) && error("No saved run to resume")
        lastdir=joinpath(dir,"segments",last(segments))
        checkpoints=sort(filter(s->startswith(s,"checkpoint_") && endswith(s,".vtkhdf"),readdir(lastdir)))
        isempty(checkpoints) && error("No checkpoint to resume")
        restart=load_checkpoint(joinpath(lastdir,last(checkpoints)),Float64,QuinticSplineMethod())
        restart.time.t >= END_TIME && error("Case already reached t=$END_TIME")
    end
    segment=joinpath(dir,"segments",@sprintf("%03d",length(segments)+1));mkdir(segment)
    write(joinpath(segment,"simulation_common_source.jl"),COMMON_SOURCE)
    write(joinpath(segment,"case_source.jl"),read(entry_file,String))
    config=Dict("case"=>name,"initial_condition"=>"VortexPasta.PredefinedCurves; uniform arclength",
        "initial_nodes"=>length.(fs0),"Gamma"=>GAMMA,"a"=>CORE,"Delta"=>DELTA,"R"=>RADIUS,
        "L"=>BOX,"mesh"=>collect(MESH),"beta"=>BETA,"quadrature_order"=>QUAD_ORDER,
        "end_time"=>END_TIME,"dt_max"=>DT_MAX,"kelvin_period_factor"=>0.25,
        "velocity_displacement_factor"=>0.15,"hmin"=>hmin,"hmax"=>hmax,
        "reconnection_distance"=>DREC,"reconnection_threads"=>1,"reconnection_max_passes"=>2,
        "reconnection_use_velocity"=>true,"reconnection_cos_max"=>0.97,
        "reconnection_decrease_length"=>true,"scheme"=>"RK4","LIA"=>false,
        "mutual_friction"=>0.0,"imposed_flow"=>false,"fold_periodic"=>false,
        "julia"=>string(VERSION),"VortexPasta"=>string(pkgversion(VortexPasta)),
        "Julia_threads"=>Threads.nthreads(),"BLAS_threads"=>BLAS.get_num_threads(),
        "manifest_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"Manifest.toml")))),
        "diagnostic_interval"=>DIAG_INTERVAL,"snapshot_interval"=>SNAPSHOT_INTERVAL,
        "checkpoint_interval"=>CHECKPOINT_INTERVAL)
    if resume
        original=TOML.parsefile(joinpath(dir,"config.toml"))
        for key in ("Gamma","a","Delta","R","L","mesh","beta","quadrature_order",
                    "initial_nodes","dt_max","kelvin_period_factor","velocity_displacement_factor",
                    "hmin","hmax","reconnection_distance","reconnection_threads",
                    "reconnection_max_passes","reconnection_use_velocity","scheme","LIA",
                    "mutual_friction","fold_periodic","manifest_sha256")
            config[key]==original[key] || error("Resume would change $key; existing parameters must be retained")
        end
    else
        save_toml(joinpath(dir,"config.toml"),config)
    end
    save_toml(joinpath(segment,"config.toml"),config)
    p=parameters()
    prob=restart===nothing ? VortexFilamentProblem(fs0,END_TIME,p) :
        VortexFilamentProblem(restart,END_TIME-restart.time.t,p)
    adapt=AdaptBasedOnVelocity(0.15min(DREC,hmin)) | MaximumTimestep(DT_MAX) |
        AdaptBasedOnSegmentLength(0.25)
    rec=ReconnectFast(DREC;nthreads=1,max_passes=2,use_velocity=true)
    started=time()
    # Keep continuous coordinates for displacement and geometric diagnostics.
    # fold_periodic controls coordinate recentering, not periodic dynamics (Ls above).
    # analyze.py wraps and splits curves into [0,2π]^3 for periodic visualization.
    iter=init(prob,RK4();dt=DT_MAX,refinement=RefineBasedOnSegmentLength(hmin,hmax),
        reconnect=rec,adaptivity=adapt,fold_periodic=false)
    open(joinpath(segment,"solver.txt"),"w") do io;show(io,iter);end
    diag=open(joinpath(segment,"diagnostics.csv"),"w")
    println(diag,"t,tstar,step,dt,loops,nodes,ds_min,ds_mean,ds_max,kappa_max,length,energy,Ix,Iy,Iz,H_over_Gamma2,cx,cy,cz,reconnections,reconnection_length_loss,refinement_length_loss,removed_loops,removed_length")
    snaps=open(joinpath(segment,"snapshots.csv"),"w");println(snaps,"t,tstar,step,file")
    events=open(joinpath(segment,"events.csv"),"w")
    println(events,"t_before,t_after,tstar_after,events_in_step,cumulative_events,loops_before,loops_after,reconnection_length_loss,refinement_length_loss,removed_loops,removed_length")
    nextdiag=iter.t;nextsnap=iter.t;nextcheck=iter.t
    lastdiag=-1;lastsnap=-1;lastcheckpoint=-1;event=false
    wall_limit=parse(Float64,get(ENV,"VORTEX_WALL_SECONDS","1200"))
    status="running";reason="";failed=nothing
    function output_diag()
        println(diag,join(diagnostic(iter),','));flush(diag);lastdiag=iter.nstep
    end
    function output_snapshot()
        file=@sprintf("curve_%08d.csv",iter.nstep)
        write_curves(joinpath(segment,file),iter.fs)
        println(snaps,join((iter.t,GAMMA*iter.t,iter.nstep,file),','));flush(snaps);lastsnap=iter.nstep
    end
    function checkpoint()
        save_checkpoint(joinpath(segment,@sprintf("checkpoint_%08d.vtkhdf",iter.nstep)),iter)
        lastcheckpoint=iter.nstep
    end
    try
        while true
            if iter.t>=nextdiag || event || iter.t>=END_TIME
                output_diag();nextdiag=iter.t+DIAG_INTERVAL
            end
            if iter.t>=nextsnap || event || iter.t>=END_TIME
                output_snapshot();nextsnap=iter.t+SNAPSHOT_INTERVAL
            end
            if iter.t>=nextcheck
                checkpoint();nextcheck=iter.t+CHECKPOINT_INTERVAL
            end
            if iter.t>=END_TIME
                status="complete";reason="Reached requested physical time";break
            elseif isempty(iter.fs)
                status="stopped";reason="No vortices remain";break
            elseif time()-started>wall_limit
                status="checkpointed";reason="Wall-time limit; --resume continues the same case";break
            elseif sum(length,iter.fs)>4000
                status="stopped";reason="Node-count guard exceeded";break
            end
            oldt=iter.t;oldrec=iter.stats.reconnection_count;oldloops=length(iter.fs)
            iter.time.dt=min(iter.dt,END_TIME-iter.t)
            step!(iter)
            event=iter.stats.reconnection_count!=oldrec
            if event
                println(events,join((oldt,iter.t,GAMMA*iter.t,iter.stats.reconnection_count-oldrec,
                    iter.stats.reconnection_count,oldloops,length(iter.fs),
                    iter.stats.reconnection_length_loss,iter.stats.refinement_length_loss,
                    iter.stats.filaments_removed_count,iter.stats.filaments_removed_length),','));flush(events)
            end
            if iter.nstep%500==0
                @printf("%s t=%.6f/%g step=%d nodes=%d loops=%d rec=%d wall=%.1fs\n",
                    name,iter.t,END_TIME,iter.nstep,sum(length,iter.fs;init=0),length(iter.fs),
                    iter.stats.reconnection_count,time()-started);flush(stdout)
            end
        end
    catch e
        status="failed";reason=sprint(showerror,e);failed=e
        open(joinpath(segment,"failure.txt"),"w") do io;showerror(io,e,catch_backtrace());end
    finally
        lastdiag!=iter.nstep && output_diag()
        lastsnap!=iter.nstep && output_snapshot()
        lastcheckpoint!=iter.nstep && checkpoint()
        close(diag);close(snaps);close(events)
        save_toml(joinpath(segment,"completion.toml"),Dict("status"=>status,"reason"=>reason,
            "final_time"=>iter.t,"requested_final_time"=>END_TIME,"elapsed_seconds"=>time()-started,
            "steps"=>iter.nstep,"nodes"=>sum(length,iter.fs;init=0),"loops"=>length(iter.fs),
            "reconnections"=>iter.stats.reconnection_count,"finished_utc"=>string(now(UTC))))
        open(joinpath(segment,"timings.txt"),"w") do io;show(io,iter.to);end
    end
    println(name,": ",status," at t=",iter.t,"; ",reason)
    failed!==nothing && throw(failed)
end
