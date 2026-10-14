include("simulate.jl")
const AUDIT_RECORDS=Vector{Vector{Float64}}()
function audit_measure(iter,stage)
 fs=iter.fs;p=iter.prob.p
 if any(f->!Filaments.check_nodes(Bool,f),fs);return;end
 cache=BiotSavart.init_cache(p)
 vs=map(similar∘nodes,fs);psi=map(similar∘nodes,fs)
 BiotSavart.compute_on_nodes!((velocity=vs,streamfunction=psi),cache,fs)
 push!(AUDIT_RECORDS,[iter.t,iter.nstep,stage,Diagnostics.kinetic_energy(fs,psi,p),Diagnostics.filament_length(fs;quad=p.quad),length(fs),sum(length,fs),iter.stats.reconnection_count])
end
include("instrument_operations.jl")
for file in ARGS
 empty!(AUDIT_RECORDS);cfg=TOML.parsefile(file);runone(cfg)
 open(joinpath(ROOT,"runs",cfg["name"],"operation_budget.csv"),"w") do io
  println(io,"t,step,stage,energy,length,loops,nodes,reconnections")
  for row in AUDIT_RECORDS;println(io,join(row,','));end
 end
end
