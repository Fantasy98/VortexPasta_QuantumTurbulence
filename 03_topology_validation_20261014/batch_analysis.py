from analyze import *
# Process only completed runs. Expensive, invalid late wide-core geometry is kept raw.
for d in sorted((ROOT/'runs').iterdir()):
 if not (d/'complete.toml').exists():continue
 if any(x in d.name for x in ['pilot','operation']):continue
 if d.name=='hopf_prepared_wide':continue
 if not (ROOT/'analysis'/d.name/'summary.json').exists():
  analyze(d);plot(d,movie=d.name in ['ring_base','hopf_base','trefoil_base','hopf_prepared_ref','trefoil_prepared_ref'])
