from analyze import *
rows=[]
for d in sorted((ROOT/'runs').iterdir()):
 if not (d/'complete.toml').exists():continue
 cfg,a,sn=loadrun(d);row={'run':d.name,'a':cfg['a'],'min_segment_over_a':float(min(a['ds_min'])/cfg['a']),'max_a_kappa':float(max(a['kappa_max'])*cfg['a']),'initial_a_kappa':float(a['kappa_max'][0]*cfg['a']),'max_nodes':int(max(a['nodes']))}
 p=ROOT/'analysis'/d.name/'geometry.json'
 if p.exists():
  g=json.loads(p.read_text());row['min_nonlocal_separation_over_a']=min(x['strand_sites'][0]['distance'] for x in g)/cfg['a']
 rows.append(row)
(ROOT/'analysis/scale_separation.json').write_text(json.dumps(rows,indent=2))
fig,axs=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
for name in ['hopf_base','trefoil_base','hopf_prepared_ref','trefoil_prepared_ref','hopf_prepared_fine','trefoil_prepared_fine']:
 d=ROOT/'runs'/name
 if not (d/'complete.toml').exists():continue
 c,a,sn=loadrun(d);axs[0].plot(a['tstar'],a['kappa_max']*c['a'],label=name)
 p=ROOT/'analysis'/name/'geometry.json'
 if p.exists():
  g=json.loads(p.read_text());axs[1].plot([x['tstar'] for x in g],[x['strand_sites'][0]['distance']/c['a'] for x in g],label=name)
axs[0].set(ylabel='a × maximum nodal curvature',yscale='log');axs[1].set(ylabel='Minimum nonlocal separation / a',yscale='log')
for ax in axs:ax.set(xlabel='t* = Γt/R²');ax.axhline(1,color='gray',ls='--');ax.grid(alpha=.2);ax.legend(fontsize=6)
fig.suptitle('Thin-core validity requires aκ ≪ 1 and nonlocal separation/a ≫ 1');fig.savefig(ROOT/'figures/scale_separation.png',dpi=160);plt.close(fig)
