include("simulate.jl")
open("analysis/evolved_accuracy.csv","w") do io
 println(io,"case,file,a,relative_velocity_L2,max_velocity_error,relative_energy_error,Hdiff")
 for name in ("hopf_base","trefoil_base","hopf_prepared_ref","trefoil_prepared_ref")
  dir=joinpath(ROOT,"runs",name);isfile(joinpath(dir,"complete.toml"))||continue
  cfg=TOML.parsefile(joinpath(dir,"config.toml"));a=cfg["a"]
  curves=sort(filter(x->startswith(x,"curve_"),readdir(dir)))
  pre=sort(filter(x->startswith(x,"pre_event_"),readdir(dir)))
  selected=isempty(pre) ? [curves[end]] : [pre[1],curves[findfirst(x->x>replace(pre[1],"pre_event_"=>"curve_"),curves)],curves[end]]
  for fname in selected
   data=readdlm(joinpath(dir,fname),',',Float64;skipstart=1)
   fs=[Filaments.init(ClosedFilament,[Vec3(row[2:4]) for row in eachrow(data[data[:,1].==i,:][1:2:end,:])],QuinticSplineMethod()) for i in unique(data[:,1])]
   p=parameters(a);pr=parameters(a;mesh=64,beta=5.,quad=12)
   it=init(VortexFilamentProblem(fs,.01,p),RK4();dt=.0001);ir=init(VortexFilamentProblem(fs,.01,pr),RK4();dt=.0001)
   dif=[norm(v-w) for (vs,ws) in zip(it.vs,ir.vs) for (v,w) in zip(vs,ws)];den=sum(sum(abs2,v) for vs in ir.vs for v in vs)
   println(io,join((name,fname,a,sqrt(sum(abs2,dif)/den),maximum(dif),abs(Diagnostics.kinetic_energy(it)/Diagnostics.kinetic_energy(ir)-1),abs(Diagnostics.helicity(it)-Diagnostics.helicity(ir))/GAMMA^2),','));flush(io)
  end
 end
end
