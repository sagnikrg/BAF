# Vision (Chapters 2, 3): simulation wishlist

Figure numbers from the full draft of 2026-11-08. All random h, the noisy circuit of the thesis, mbl_dtc code.

### V1: refined epsilon windows at L = 8, 10, 12 (kind 1, option B)
- **Status:** queued on BAF 2026-11-10 (5412 jobs: L=8,10 all 7 theta, 4000 real.; L=12 theta 0 x1000, 0.1/0.2/0.3 x500; L=12 theta 0.05/0.15/0.25 later). Test passed (E:/mbl_dtc/baf_test)
- **Where:** cluster (BAF)
- **Parameters:** new g points at the midpoints of the old 0.025 grid, in windows of +-0.025 around every LSR and S/L crossing per side: theta 0 (18 points), 0.05 (11), 0.1 (19), 0.15 (9), 0.2 (8), 0.25 (4, plus edge points g=0.0125, 0.9875), 0.3 (edge points only). Full lists in memory `baf-eps-plan-2026-11-07`. Window points at L=12 (and L=10 where needed) get a new uniform, larger realisation count, kept as their own dataset (no mixing of counts within a curve). DTC-side windows to be widened to cover the median band gap and eta crossings (g 0.75-0.97).
- **Cost:** about 14k core-h at L=12 (theta 0 alone ~7k), ~1.4k at L=10, ~0.1k at L=8; many short jobs (L=12 ~350-1400 s per realisation and point).
- **For:** the crossings and projected boundaries of Fig. 2.4 `fig:lsr_phase_portrait`, Fig. 3.5 `fig:phase_portraits` and Fig. 3.24 `fig:phase_portraits_gap`; the crossing errors are set by realisations (sigma_c 0.013-0.05), so this sharpens all three.
- **Notes:** the closing points at theta >= 0.1 rest on ~500 realisations at L=12 today; those rows first.

### V1.2: more L=12 realisations, window points and the old grid
- **Status:** scripted 2026-11-10 (Tony), mbl_dtc/jobs_v1.2.txt, 5500 jobs; local test passed (julia 1.9.4, both sets)
- **Parameters:** L=12, window points (set win): +1000 at theta 0 (jobs 1001-2000), +500 at 0.1, 0.2, 0.3 (jobs 501-1000) -> 2000 / 1000 per point. Old 0.025 grid (set grid, 41 points, eps/grid.txt): +500 at theta 0.1, 0.15, 0.2, 0.25, 0.3, 0.35 (jobs 501-1000) -> ~1000 each; theta 0 stays at 4066, 0.05 at 1996.
- **Cost:** ~30k core-h: windows ~6k (theta 0 ~3.5k), grid ~24k (3000 jobs x 41 pts x ~700 s, ~8 h per job, worst ~16 h)

### V2: L = 14, 16 with POLFED, two antipodal bands (kind 2)
- **Status:** approved 2026-11-10 (Tony): old theta 0, 0.1, 0.2, 0.3 first, ~20 eps per theta, L=14 x500, L=16 x100 realisations; scripting in BAF/V2
- **Where:** cluster (BAF)
- **Parameters:** L=14, 16 (L=18 postponed), one non-uniform epsilon list per theta, theta 0-0.3; fewer realisations than L=12; two bands per realisation (around E and E+pi); same J, h per realisation across epsilon from a disorder table made up front; jobs split over epsilon.
- **Cost:** ~14k core-h (L=14 ~4.7k at ~7-8 min per realisation and eps, L=16 ~9.3k at ~70 min, BAF = 2x laptop); L=14 jobs = 1 realisation x all eps of a theta (~2.5 h), L=16 jobs = 1 realisation x ~10 eps (~12 h).
- **For:** the band gap and eta crossings (Fig. 3.23 `fig:band_gap_median`, Fig. 3.20 `fig:psd_scalar`) and the time crystal boundaries of Fig. 3.24 `fig:phase_portraits_gap`, against Fig. 2.4 `fig:lsr_phase_portrait`. This decides whether the DTC boundary is distinct from the MBL one (the backbone of the thesis and of the flagship paper).
- **Notes:** the longest share of the cluster; only quantities the thesis plots.
