Three VortexPasta base cases: 64^3 periodic mesh, physical time 0 <= t <= 10
=======================================================================

All three local runs reached t=10. See RUN_SUMMARY.txt for actual execution
times and outcomes, and validation/completed_runs.json for the final checks.

This directory contains exactly three physical configurations:
  single_ring.jl   One circular ring, initially 64 filament nodes.
  hopf_link.jl     Two perpendicular linked circles, initially 64 nodes each.
  trefoil.jl       One TrefoilKnot, initially 128 nodes.

Each case constructs its initial geometry with VortexPasta.PredefinedCurves.
No GPE-prepared centerlines or CSV inputs are imported. The common Julia file
is a shared numerical/output implementation, not an additional case. There
are no fine-grid, half-timestep, wide-core, or alternative-geometry variants.

Run from this directory:
  cd /home/yuninw/codes/simulations/vortexPasta/04_helicity_64_three_cases
  ./run.sh single_ring
  ./run.sh hopf_link
  ./run.sh trefoil

The three commands above are sequential. Each defaults to two Julia threads
and one BLAS/OpenMP thread. At most two simultaneous runs are recommended on
this host. When only one case is active, four Julia threads can be selected:
  VORTEX_THREADS=4 ./run.sh trefoil --resume
Each execution segment records its actual thread count. This changes CPU
parallelism, not physical parameters or the numerical resolution.
Direct Julia commands also work:
  JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 julia --project=. single_ring.jl
  JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 julia --project=. hopf_link.jl
  JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 julia --project=. trefoil.jl

Results are separate:
  results/single_ring/
  results/hopf_link/
  results/trefoil/

A fresh launch refuses to overwrite existing results. It writes config.toml
and a segments/001/ directory containing diagnostics.csv, snapshots.csv,
events.csv, curve_*.csv, checkpoint_*.vtkhdf, source copies, timings.txt,
and completion.toml. A status of complete is written only after reaching
t=10. A bounded stop is explicitly recorded as checkpointed or stopped.

The default wall-time guard is 1200 seconds per invocation. If it is reached,
continue the same case, preserving all earlier data, with:
  ./run.sh hopf_link --resume
The continuation loads the latest full solver checkpoint and writes a new
segments/002/ directory; solver time, counters, and physical parameters are
retained. The guard can be changed explicitly through VORTEX_WALL_SECONDS.
An already completed case cannot be resumed. --check validates initialization
without running dynamics or creating a result case.

Analysis (Python 3.11+, NumPy, SciPy, Matplotlib, ffmpeg):
  OPENBLAS_NUM_THREADS=1 python3 analyze.py single_ring
  OPENBLAS_NUM_THREADS=1 python3 analyze.py hopf_link
  OPENBLAS_NUM_THREADS=1 python3 analyze.py trefoil
or:
  OPENBLAS_NUM_THREADS=1 python3 analyze.py all

The script selects case-specific plots and writes results/CASE/analysis/:
  diagnostics.png / .pdf   Helicity, length, energy, impulse, loops, nodes,
                          segment lengths, curvature, numerical length losses.
  case_analysis.png/.pdf   Ring displacement/speed, or linking/writhe/events.
  snapshots.png / .pdf     Four physical-time snapshots.
  evolution.mp4           Fixed laboratory frame, equal physical scales,
                          physical and dimensionless timestamps.
  geometry.json           Approximate Gauss linking/writhe histories.
  summary.json            Actual endpoint, counts, checks, and interpretation.

Analysis refuses to overwrite an existing analysis folder. To reanalyze with
new plotting options or code, select a fresh output directory:
  python3 analyze.py single_ring --output results/single_ring/analysis_revision2
Use --no-video to generate plots only. Videos have approximately 20 seconds
of playback; frame selection is approximately uniform in physical time.
The ring uses unwrapped laboratory z and retains its initial position as a
gray reference. Its movement is also shown on a separate displacement panel.
Movies do not subtract the moving centroid.

Parameters and conventions
--------------------------
These are the previous study's slender-core base parameters, with only the
Fourier mesh changed from 48^3 to 64^3 and the requested endpoint set to t=10:
  Periodic box side L = 2*pi; R = 1.
  Gamma = 2.0826013772617338; a = 0.001; Delta = 0.615.
  Full periodic Ewald Biot-Savart, CPU, RK4, no imposed flow, no mutual friction.
  Ewald beta = 4.5; quadrature order = 10 (normal and near-singular).
  NUFFT oversampling = 1.5, half support = 6; LIA segment fraction = 0.2.
  Maximum timestep = 0.001, Kelvin-period factor = 0.25, displacement factor = 0.15.
  Quintic splines, approximate uniform arclength initialization, adaptive nodes.
  Refinement bounds = 0.65 and 1.4 times the smallest initial segment length.
  ReconnectFast distance = 0.08, velocity check on, one reconnection thread,
  at most two scans per step, package defaults cos_max=0.97 and decrease_length=true.
  Diagnostics every ~0.02, curve snapshots every ~0.05, checkpoints every ~1.
  Additional diagnostics and snapshots are saved after reconnection steps.

64^3 is the Fourier/Ewald mesh, not 64^3 filament nodes or a GPE density grid.
The actual evolving filament node count and segment lengths are in diagnostics.
The full t=10 interval corresponds to t*=Gamma*t/R^2=20.8260137726.

The ring starts at (pi,pi,pi), oriented counterclockwise in the xy plane.
The Hopf centers are (pi-1/2,pi,pi) and (pi+1/2,pi,pi); the second Ring is
rotated +pi/2 around x, giving oriented linking number -1.
TrefoilKnot is scaled by R/3 and translated to (pi,pi,pi), orientation +1.
Arclength resampling changes the parameter-to-node mapping, not the curve.

Native helicity is H=Gamma*integral(v_BS.ds), reported as H/Gamma^2. The
geometric crosscheck is sum(Wr)+2sum(Lk), without ribbon twist. Its midpoint
Gauss quadrature on resampled exported curves is approximate, especially
near reconnection. Noncontractible periodic loops are explicitly excluded
from that Euclidean Gauss diagnostic. It uses a common continuous lift and
does not sum periodic-image Gauss terms; its interpretation is limited for
extended curves. Per-loop coordinate spans are saved alongside it. Native
velocity-based helicity uses the full periodic solver. No density, sound, or GPE spectrum is
computed. Energy is E_total/(rho*V), impulse is I/(rho*Gamma).
Reconnection, refinement and unresolved-loop removal are numerical operations;
their recorded length changes are not identified as physical dissipation.

These three runs do not supply a new convergence study or a quantitative
comparison to published paper curves. No GPE solver data are plotted here.
Project.toml and Manifest.toml pin the inherited Julia environment. The tested
runtime is Julia 1.13.1 / VortexPasta 0.34.19; Python versions are recorded in
requirements-analysis.txt. GLMakie remains in the inherited environment but
is not loaded by these scripts. No remote jobs or paid services are used.
