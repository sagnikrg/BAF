# Ultron (Chapter 8, and the TEV benchmark and POLFED sections of Chapter 2): simulation wishlist

Figure numbers from the full draft of 2026-11-08 (main.aux, 14:15). Nothing quick and nothing urgent: both parts are written on data that is already there. The open entries would turn statements of the text into measured numbers.

### U1: power basis at L = 18 with the matrix product state evolution, to the end
- **Status:** planned, 2026-11-08
- **Where:** cluster (BAF), one chain of jobs
- **Parameters:** L=18, theta=0, g=1, one realisation, cutoff 1e-14, n_cv=64, k*n_cv = 2d (N = 524289 steps), script `polfed/code/polfed.large.jl` (checkpoint every 20000 steps, resumes on rerun)
- **Cost:** the first 2e4 steps took 1.25 s per step at chi ~ 490 and growing, so >= 180 core-h, i.e. >= 4 jobs of 50 h resumed from the checkpoint; memory < 2 GiB (Krylov part ~1.5 GiB)
- **For:** Fig. 8.10 `fig:polfed_scaling`, Tab. 8.3 `tab:polfed_scaling` and Sec. `sec:polfed_reach`: today the L=18 point is "stopped after 2e4 steps, more than a week"; the run would give the measured time, bond dimension and accuracy instead.
- **Notes:** low priority, the conclusion of the section does not change. The checkpoint of the stopped laptop run (`polfed.large.L18.checkpoint.jls`, 2e4 steps) can be the start. Needs the cluster test of U0.

### U2: pi pairs at L = 16 from exact overlaps
- **Status:** planned, 2026-11-08
- **Where:** laptop (job queue) or cluster, single thread
- **Parameters:** L=16, theta=0, g=0.95 and 0.8, 4 realisations each, dense (exact) evolution of N = 3d = 196608 steps, k=d/16, n_cv=48, targets 0 and pi (`polfed/code/polfed.pairs.jl` with L=16 added)
- **Cost:** a dense step at L=16 is ~0.04 s, so ~2-2.5 h per realisation, 8 jobs, ~20 core-h; memory < 1 GiB
- **For:** Fig. 8.6 `fig:polfed_pairs` (now L=8-14): whether the splitting at L=16 is still resolved, or sits at the floor of ~1e-10 set by the accuracy of the eigenvalues (the median at g=0.95 is ~3e-10 at L=14 and falls with L).
- **Notes:** also the direct check of the resolution issue for Vision's V2 (see below), so if V2 goes ahead this is worth doing first.

### U3: threads for long runs at L = 24 (adaptive time evolution)
- **Status:** planned, optional, 2026-11-08
- **Where:** laptop in the office (parallel, so not at home) or cluster with 4 cores
- **Parameters:** L=24, theta=0, g=0.9 and 0.95, 500 Floquet steps, AdaptiveTEV with 1 and 4 threads, 2 realisations
- **Cost:** ~10-30 core-h, longest job ~10 h
- **For:** Fig. 2.26 `fig:adaptive_tev_threads`: the text says the cores only start to matter for longer runs or larger L; this would show where.
- **Notes:** nice to have only.

### U0: test of the polfed module on the cluster environment
- **Status:** planned, 2026-11-08
- **Where:** cluster (BAF)
- **Parameters:** L=8, 10, a few points, normal POLFED and the power basis against full diagonalisation of the same realisation
- **Cost:** < 1 core-h, minutes per job
- **For:** nothing in the thesis directly; U1 and Vision's V2 depend on it.
- **Notes:** the polfed module was written and tested with julia 1.11/1.12, ITensors 0.9.11, ITensorMPS 0.4.1; the cluster depot (`.julia/julia-1.9.4-08-08-24.tar.gz`) has julia 1.9.4, ITensors 0.6.16 and Arpack. Whether it runs there unchanged is not known.

## Done (final)

### F1: Chapter 2 time evolution benchmarks
- **Status:** final (`half_chain_ee/code/{equivalence,benchmark,scaling,scaling.long,scaling.map,accuracy,parallel,parallel.chi,dynamic,adaptive}.tev.hdf5`), 2026-11-03 to 2026-11-05
- **Where:** laptop, 7 single threaded workers
- **For:** Fig. 2.14-2.26 (`fig:equivalence_tev` to `fig:adaptive_tev_threads`), Tab. 2.2 `tab:benchmark_tev`, Tab. 2.3 `tab:scaling_tev`

### F2: Chapter 2 POLFED (dense, normal)
- **Status:** final (`polfed/code/polfed.{validation,parameters,scaling,entanglement}.hdf5`), 2026-11-03/04
- **Where:** laptop
- **For:** Fig. 2.27-2.31 (`fig:polfed_filter` to `fig:polfed_entanglement`)

### F3: Chapter 8 power basis
- **Status:** final (`polfed/code/polfed.{example,krylov,krylov.overlaps,krylov.breakdown,krylov.gpu,krylov.gpu.large,precision,pairs,polynomial,polynomial.window,validation,accuracy,scaling}.hdf5`), 2026-11-04/05
- **Where:** laptop (CPU, and the GPU for `tab:polfed_gpu`)
- **For:** Fig. 8.1-8.10, Tab. 8.1-8.3

## Notes on Vision's V2 (L = 14, 16 with POLFED, two antipodal bands)

It runs on my code, so here is what I see. None of it blocks it, but 2 and 4 should be settled before scripting.

1. **Cost.** Normal POLFED (dense, ARPACK) took 3 min at L=14 and 30 min at L=16 per target on the laptop, n_ev=20, memory < 1 GiB (Fig. 2.30 `fig:polfed_scaling_normal`), growing ~8-11x per two sites (L=18: ~5 h and ~4 GiB per target). Two bands are two runs, so ~1 h per realisation and point at L=16. This was measured at theta=0, g=1; more eigenvalues per band need a larger n_cv, and the restarts may take longer elsewhere in the phase diagram. A 36 h job is then ~30 realisations-and-points at L=16. The power basis is not the better tool here: it is only cheap near theta=0, g=1, and with exact overlaps (N=3d steps, ~2 h at L=16) it costs about as much as the two ARPACK runs.
2. **Resolution of the splitting.** The eigenvalues are accurate to ~1e-9-1e-10. If the band gap is the pi-pair splitting, deep in the time crystal it already sits at that floor at L=14 (median ~3e-10 at g=0.95) and below it at L=16. Near the boundary the splittings are much larger and the crossings are not affected, but deep points should be marked as "below resolution" rather than plotted as values. U2 measures exactly this.
3. **Pairing across two runs.** The two runs give two windows; eigenvalues near the edge of one window have their partner outside the other, so the edges have to be dropped (Fig. 8.6 used a tolerance of a tenth of the level spacing).
4. **Disorder table and environment.** The ITensor circuit `circuit_dtc(sites, theta, epsilon)` draws J and h itself; a method taking given J and h (as the matrix version already has) is a small addition to the header, mine to do under the lock when V2 is scripted. The noise is drawn when the circuit is built and stays fixed during the run, consistent with not storing it. And the cluster environment is untested for this code (U0).
