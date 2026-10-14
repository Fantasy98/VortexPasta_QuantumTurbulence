using VortexPasta
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using VortexPasta.PredefinedCurves
using LinearAlgebra, Statistics, DelimitedFiles, TOML, SHA, Printf
LinearAlgebra.BLAS.set_num_threads(1)
const ROOT=@__DIR__
const GAMMA=4π*(2*0.1171875/sqrt(2))
const EVENT_ROWS=Vector{Vector{Float64}}()
const EVENT_TIME=Ref(0.0)
function record_event(nodes,is,js,Ls,lb,la)
    x=nodes[is[1]]; pts=[x+[mod(v[k]-x[k]+Ls[k]/2,Ls[k])-Ls[k]/2 for k=1:3] for v in (nodes[is[2]],nodes[js[1]],nodes[js[2]])]
    c=(x+sum(pts))/4
    push!(EVENT_ROWS,[EVENT_TIME[],c...,lb,la,is...,js...])
end
# Exact installed serial reconnection routine, with a read-only event logging hook.
include("instrument_reconnections.jl")
function parameters(a; mesh=48,beta=4.5,quad=10,periodic=true)
    common=(Γ=GAMMA,a=a,Δ=.615,quadrature=GaussLegendre(quad),quadrature_near_singularity=GaussLegendre(quad),lia_segment_fraction=.2)
    if !periodic
        return ParamsBiotSavart(;common...,Ls=Infinity(),α=Zero(),backend_short=NaiveShortRangeBackend())
    end
    alpha=(mesh-1)/(4beta);rcut=beta/alpha
    ParamsBiotSavart(;common...,Ls=(2π,2π,2π),Ns=(mesh,mesh,mesh),α=alpha,rcut,
        backend_short=CellListsBackend(2),backend_long=NonuniformFFTsBackend(CPU();σ=1.5,m=HalfSupport(6)))
end
function getfs(input)
    data=readdlm(input,',',Float64;skipstart=1)
    [Filaments.init(ClosedFilament, [Vec3(row[2:4]) for row in eachrow(data) if row[1]==i],QuinticSplineMethod()) for i in unique(data[:,1])]
end
function writecurves(file,fs; refine=1)
    open(file,"w") do io
        println(io,"loop,x,y,z")
        for (j,f) in enumerate(fs), i in eachindex(f), q=0:refine-1
            v=f(i,q/refine); println(io,join((j,v...),','))
        end
    end
end
function features(iter)
    fs=iter.fs;p=iter.prob.p
    ds=[norm(f[i+1]-f[i]) for f in fs for i in eachindex(f)]
    ks=[f[i,CurvatureScalar()] for f in fs for i in eachindex(f)]
    pos=[x for f in fs for x in nodes(f)]; c=mean(pos)
    mom=Diagnostics.vortex_impulse(iter)
    [iter.t,GAMMA*iter.t,iter.nstep,iter.dt,length(fs),sum(length,fs),minimum(ds),mean(ds),maximum(ds),maximum(ks),
      Diagnostics.filament_length(iter),Diagnostics.kinetic_energy(iter),mom...,Diagnostics.helicity(iter)/GAMMA^2,c...,
      iter.stats.reconnection_count,iter.stats.reconnection_length_loss,iter.stats.refinement_length_loss,
      iter.stats.filaments_removed_count,iter.stats.filaments_removed_length]
end
function runone(cfg)
    name=cfg["name"];dir=joinpath(ROOT,"runs",name);ispath(dir)&&error("Refusing to overwrite $dir");mkdir(dir)
    cp(@__FILE__,joinpath(dir,"simulate_source.jl"));cp(joinpath(ROOT,"instrument_reconnections.jl"),joinpath(dir,"instrument_source.jl"))
    input=joinpath(ROOT,"inputs",cfg["input"]);cp(input,joinpath(dir,"initial.csv"));cfg["input_sha256"]=bytes2hex(sha256(read(input)))
    a=get(cfg,"a",.001);d=get(cfg,"drec",.08);dt=get(cfg,"dt",.001);tf=cfg["end_time"]
    fs=getfs(input);h=minimum(norm(f[i+1]-f[i]) for f in fs for i in eachindex(f))
    hmin=.65h;hmax=1.4h
    p=parameters(a;mesh=get(cfg,"mesh",48),beta=get(cfg,"beta",4.5),quad=get(cfg,"quad",10))
    cfg=merge(cfg,Dict("Gamma"=>GAMMA,"Delta"=>.615,"R"=>1.,"L"=>2π,"a"=>a,"hmin"=>hmin,"hmax"=>hmax,
      "julia"=>string(VERSION),"VortexPasta"=>string(pkgversion(VortexPasta)),"threads"=>Threads.nthreads(),"mesh"=>get(cfg,"mesh",48),"beta"=>get(cfg,"beta",4.5),"quad"=>get(cfg,"quad",10),"mutual_friction"=>0.,"LIA"=>false,
      "reconnection_threads"=>1,"reconnection_max_passes"=>2,"reconnection_cos_max"=>.97,"decrease_length"=>true,"use_velocity"=>true))
    open(joinpath(dir,"config.toml"),"w") do io;TOML.print(io,cfg);end
    prob=VortexFilamentProblem(fs,tf,p)
    adapt=AdaptBasedOnVelocity(get(cfg,"velocity_factor",.15)*min(d,hmin))|MaximumTimestep(dt)
    # The thin-core Kelvin-period formula is invalid for the wide-core exploratory case.
    if a<.01; adapt=adapt|AdaptBasedOnSegmentLength(get(cfg,"dtfactor",.25));end
    rec=ReconnectFast(d;nthreads=1,max_passes=2,use_velocity=true)
    empty!(EVENT_ROWS)
    started=time();iter=init(prob,RK4();dt,reconnect=rec,refinement=RefineBasedOnSegmentLength(hmin,hmax),adaptivity=adapt,fold_periodic=false)
    open(joinpath(dir,"solver.txt"),"w") do io;show(io,iter);end
    headers="t,tstar,step,dt,loops,nodes,ds_min,ds_mean,ds_max,kappa_max,length,energy,Ix,Iy,Iz,H_over_Gamma2,cx,cy,cz,reconnections,reconnection_length_loss,refinement_length_loss,removed_loops,removed_length"
    diag=open(joinpath(dir,"diagnostics.csv"),"w");println(diag,headers)
    snaps=open(joinpath(dir,"snapshots.csv"),"w");println(snaps,"t,tstar,step,file")
    eventfile=open(joinpath(dir,"events.csv"),"w");println(eventfile,"t_before,t_after,tstar_after,x,y,z,length_before,length_after,i1,i2,j1,j2")
    nextdiag=0.;nextsnap=0.;nextcheck=.5;lastrec=0;lastout=-1;status="complete";stopreason="target reached"
    function output_snapshot(iter,suffix="")
        fname=@sprintf("curve_%07d%s.csv",iter.nstep,suffix);writecurves(joinpath(dir,fname),iter.fs;refine=2)
        println(snaps,join((iter.t,GAMMA*iter.t,iter.nstep,fname),','));flush(snaps)
    end
    save_checkpoint(joinpath(dir,"checkpoint_initial.vtkhdf"),iter)
    while true
        if iter.t+1e-12>=nextdiag || iter.stats.reconnection_count!=lastrec || iter.t>=tf
            println(diag,join(features(iter),','));flush(diag);nextdiag=iter.t+get(cfg,"diagnostic_interval",.005)
        end
        if iter.t+1e-12>=nextsnap || iter.stats.reconnection_count!=lastrec || iter.t>=tf
            output_snapshot(iter);nextsnap=iter.t+get(cfg,"snapshot_interval",.02);lastout=iter.nstep
        end
        if iter.t>=nextcheck
            save_checkpoint(joinpath(dir,@sprintf("checkpoint_%07d.vtkhdf",iter.nstep)),iter);nextcheck=iter.t+.5
        end
        lastrec=iter.stats.reconnection_count
        if iter.t>=tf;break;end
        if time()-started>get(cfg,"wall_limit",900.) || sum(length,iter.fs)>4000 || isempty(iter.fs)
            status="bounded_stop";stopreason="wall/node/empty-filament guard";break
        end
        before=[copy(f) for f in iter.fs];oldt=iter.t;oldstep=iter.nstep
        EVENT_TIME[]=oldt
        step!(iter)
        if !isempty(EVENT_ROWS)
            writecurves(joinpath(dir,@sprintf("pre_event_%07d.csv",oldstep)),before;refine=2)
            for row in EVENT_ROWS
                println(eventfile,join((row[1],iter.t,GAMMA*iter.t,row[2:end]...),','))
            end
            flush(eventfile);empty!(EVENT_ROWS)
        end
        if iter.nstep%500==0
            @printf("%s step=%d t=%.5f/%g nodes=%d loops=%d rec=%d wall=%.1fs\n",name,iter.nstep,iter.t,tf,sum(length,iter.fs),length(iter.fs),iter.stats.reconnection_count,time()-started);flush(stdout)
        end
    end
    lastout!=iter.nstep&&output_snapshot(iter)
    save_checkpoint(joinpath(dir,"checkpoint_final.vtkhdf"),iter)
    close(diag);close(snaps);close(eventfile)
    summary=Dict("status"=>status,"reason"=>stopreason,"elapsed_seconds"=>time()-started,"final_time"=>iter.t,"final_tstar"=>GAMMA*iter.t,"steps"=>iter.nstep,"nodes"=>sum(length,iter.fs),"loops"=>length(iter.fs),"reconnections"=>iter.stats.reconnection_count)
    open(joinpath(dir,"complete.toml"),"w") do io;TOML.print(io,summary);end
    open(joinpath(dir,"timings.txt"),"w") do io;show(io,iter.to);end
    println(name," ",summary);flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__
    for file in ARGS
        cfg=TOML.parsefile(file)
        try;runone(cfg)
        catch e
            dir=joinpath(ROOT,"runs",cfg["name"]);mkpath(dir)
            open(joinpath(dir,"failed.txt"),"w") do io;showerror(io,e,catch_backtrace());end
            showerror(stderr,e,catch_backtrace());println(stderr);rethrow()
        end
    end
end
