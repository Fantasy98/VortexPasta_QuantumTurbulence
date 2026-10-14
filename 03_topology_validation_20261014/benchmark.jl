include("simulate.jl")
open(joinpath(ROOT,"analysis/ring_benchmark.csv"),"w") do io
 println(io,"a,n,periodic,mesh,beta,U,Upred,E,Epred,H_over_Gamma2,Ix,Iy,Iz")
 for a in (.001,.1171875), n in (64,128), periodic in (false,true)
  fs=getfs(joinpath(ROOT,"inputs/ring_$(n).csv"));p=parameters(a;periodic)
  iter=init(VortexFilamentProblem(fs,.01,p),RK4();dt=.0001)
  u=mean(x[3] for v in iter.vs for x in v)
  pred=GAMMA/(4π)*(log(8/a)-.615);epred=GAMMA^2/2*(log(8/a)-1.615)
  println(io,join((a,n,periodic,32,3.5,u,pred,Diagnostics.kinetic_energy(iter),epred,Diagnostics.helicity(iter)/GAMMA^2,Diagnostics.vortex_impulse(iter)...),','));flush(io)
 end
end
# Independently tighten Ewald splitting, mesh, quadrature and NUFFT support.
open(joinpath(ROOT,"analysis/biot_savart_accuracy.csv"),"w") do io
 println(io,"input,a,relative_velocity_L2,max_velocity_error,relative_energy_error,absolute_H_error")
 for input in ("ring_128.csv","hopf_128.csv","trefoil_256.csv","hopf_prepared_128.csv","trefoil_prepared_256.csv"), a in (.001,.1171875)
  fs=getfs(joinpath(ROOT,"inputs",input));p=parameters(a);pref=parameters(a;mesh=64,beta=4.5,quad=6)
  it=init(VortexFilamentProblem(fs,.01,p),RK4();dt=.0001)
  ir=init(VortexFilamentProblem(fs,.01,pref),RK4();dt=.0001)
  dif=[norm(v-w) for (vs,ws) in zip(it.vs,ir.vs) for (v,w) in zip(vs,ws)]
  den=sum(sum(abs2,v) for vs in ir.vs for v in vs)
  println(io,join((input,a,sqrt(sum(abs2,dif)/den),maximum(dif),abs(Diagnostics.kinetic_energy(it)/Diagnostics.kinetic_energy(ir)-1),abs(Diagnostics.helicity(it)-Diagnostics.helicity(ir))/GAMMA^2),','));flush(io)
 end
end
S=define_curve(TrefoilKnot();scale=1/3,translate=π,orientation=1)
err=maximum(norm(S(t)-Vec3(π+(sin(2π*t)+2sin(4π*t))/3,π+(cos(2π*t)-2cos(4π*t))/3,π-sin(6π*t)/3)) for t in range(0,1,length=1001))
@assert err<1e-14
println("Trefoil definition maximum coordinate error=",err)
