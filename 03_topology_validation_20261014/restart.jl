# Usage: julia --project=. restart.jl runs/CASE/checkpoint_final.vtkhdf inputs/NEW_RESTART.toml
# The new config must explicitly preserve a, Delta, Gamma, mesh, and refinement scales.
# Export the checkpoint into a new input and start from its physical time.
include("simulate.jl")
file,config=ARGS
cfg=TOML.parsefile(config)
cp=load_checkpoint(file,Float64,QuinticSplineMethod())
newinput=cfg["name"]*"_restart.csv";dest=joinpath(ROOT,"inputs",newinput)
ispath(dest)&&error("Refusing existing restart input")
writecurves(dest,cp.fs)
cfg["input"]=newinput;cfg["start_time"]=cp.time.t
cfg["restart_source"]=abspath(file);cfg["restart_source_sha256"]=bytes2hex(sha256(read(file)))
# This restart driver restores the geometry and physical clock. Cumulative counters restart
# at zero; the exact-state load_checkpoint API is separately tested in benchmark_final.jl.
runone(cfg)
