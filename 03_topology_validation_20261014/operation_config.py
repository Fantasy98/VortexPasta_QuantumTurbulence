from pathlib import Path
import csv,tomllib,numpy as np
r=Path(__file__).resolve().parent
for case in ['hopf','trefoil']:
 d=r/'runs'/f'{case}_base'
 if not (d/'complete.toml').exists():continue
 cfg=tomllib.loads((d/'config.toml').read_text());ev=list(csv.DictReader((d/'events.csv').open()))
 if not ev:continue
 first=float(ev[0]['t_before']);last=float(ev[min(2,len(ev)-1)]['t_after'])
 snaps=sorted(d.glob('pre_event_*.csv'))
 p=snaps[0];a=np.loadtxt(p,delimiter=',',skiprows=1);parts=[]
 # Saved curves have 2 interpolation points per node; restore original node locations.
 for loop in np.unique(a[:,0]):parts.append(a[a[:,0]==loop][::2])
 name=case+'_operation_audit';input=name+'.csv';np.savetxt(r/'inputs'/input,np.concatenate(parts),delimiter=',',header='loop,x,y,z',comments='')
 out={k:cfg[k] for k in ['a','drec','dt','hmin','hmax','mesh','beta','quad']}
 out.update(name=name,input=input,start_time=first,end_time=last+.002,wall_limit=120.,diagnostic_interval=.0001,snapshot_interval=.001)
 (r/'inputs'/f'{name}.toml').write_text('\n'.join(f'{k} = '+('"'+v+'"' if isinstance(v,str) else str(v)) for k,v in out.items())+'\n')
