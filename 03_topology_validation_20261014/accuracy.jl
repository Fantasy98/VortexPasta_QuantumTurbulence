include("simulate.jl")
open("analysis/accuracy_components.csv","w") do io
 println(io,"curve,n,a,mesh,beta,quad,U,relative_velocity_L2,relative_energy_error,Hdiff")
 for (input,n) in (("ring",128),("hopf",128),("trefoil",256),("trefoil_prepared",256))
  fs=getfs("inputs/$(input)_$(n).csv")
  a=.001;pref=parameters(a;mesh=64,beta=5.,quad=12);ir=init(VortexFilamentProblem(fs,.01,pref),RK4();dt=.0001)
  for (mesh,beta,quad) in ((32,3.5,4),(32,3.5,8),(48,4.5,8),(48,4.5,10),(64,5.,12))
   p=parameters(a;mesh,beta,quad);it=init(VortexFilamentProblem(fs,.01,p),RK4();dt=.0001)
   diff=sum(sum(abs2,v-w) for (vs,ws) in zip(it.vs,ir.vs) for (v,w) in zip(vs,ws));den=sum(sum(abs2,v) for vs in ir.vs for v in vs)
   u=mean(x[3] for v in it.vs for x in v)
   println(io,join((input,n,a,mesh,beta,quad,u,sqrt(diff/den),abs(Diagnostics.kinetic_energy(it)/Diagnostics.kinetic_energy(ir)-1),abs(Diagnostics.helicity(it)-Diagnostics.helicity(ir))/GAMMA^2),','));flush(io)
  end
 end
end
