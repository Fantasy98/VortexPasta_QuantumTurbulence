include("simulate.jl")
open("analysis/ring_benchmark_final.csv","w") do io
 println(io,"a,n,periodic,U,Upred,E_total_per_density,Epred,H_over_Gamma2,impulse_over_rhoGamma_z")
 for a in (.001,.1171875), n in (64,128), periodic in (false,true)
  fs=getfs("inputs/ring_$(n).csv");p=parameters(a;periodic)
  it=init(VortexFilamentProblem(fs,.01,p),RK4();dt=.0001)
  println(io,join((a,n,periodic,mean(v[3] for vs in it.vs for v in vs),GAMMA/(4π)*(log(8/a)-.615),Diagnostics.kinetic_energy(it)*(periodic ? (2π)^3 : 1),GAMMA^2/2*(log(8/a)-1.615),Diagnostics.helicity(it)/GAMMA^2,Diagnostics.vortex_impulse(it)[3]),','));flush(io)
 end
end
# Checkpoint reconstruction and actual continuation with all parameters explicit.
fs=getfs("inputs/ring_64.csv");p=parameters(.001)
it=init(VortexFilamentProblem(fs,.02,p),RK4();dt=.0002,fold_periodic=false)
for _=1:10;step!(it);end
save_checkpoint("analysis/restart_control.vtkhdf",it)
cp=load_checkpoint("analysis/restart_control.vtkhdf",Float64,QuinticSplineMethod())
ri=init(VortexFilamentProblem(cp,.02,p),RK4();dt=.0002,fold_periodic=false)
for _=1:10;step!(it);step!(ri);end
err=maximum(norm(x-y) for (f,g) in zip(it.fs,ri.fs) for (x,y) in zip(f,g))
@assert err<1e-11
open("analysis/restart_check.toml","w") do io;TOML.print(io,Dict("max_node_difference"=>err,"time_difference"=>abs(it.t-ri.t),"steps_difference"=>it.nstep-ri.nstep));end
println("Restart continuation difference: ",err)
