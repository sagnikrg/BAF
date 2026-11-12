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
- **Status:** scripted 2026-11-10 (Tony), mbl_dtc/jobs_v1.2.txt, 4000 jobs; local test passed (julia 1.9.4, both sets)
- **Parameters:** L=12, window points (set win): +1000 at theta 0 (jobs 1001-2000), +500 at 0.1, 0.2, 0.3 (jobs 501-1000) -> 2000 / 1000 per point. Old 0.025 grid (set grid, 41 points, eps/grid.txt): +500 at theta 0.1, 0.2, 0.3 (jobs 501-1000) -> ~1000 each, matching their window points (no apples and oranges); theta 0 stays at 4066 on the grid, 2000 on the windows.
- **Cost:** ~18k core-h: windows ~6k (theta 0 ~3.5k), grid ~12k (1500 jobs x 41 pts x ~700 s, ~8 h per job, worst ~16 h)

### V2: L = 14, 16 with POLFED, two antipodal bands (kind 2)
- **Status:** approved 2026-11-10 (Tony): old theta 0, 0.1, 0.2, 0.3 first, ~20 eps per theta, L=14 x500, L=16 x100 realisations. Scripted in BAF/mbl_dtc_polfed (polfed_bands.jl, two bands at phi=0, pi, nev=64; checked against full diagonalisation at L=8, exact to 1e-14, inner half of each band complete). Waits for: the tarball (julia 1.11.7, Tony), the eps lists from the full L=8-12 data, then jobs_test.txt on BAF
- **Where:** cluster (BAF)
- **Parameters:** L=14, 16 (L=18 postponed), one non-uniform epsilon list per theta, theta 0-0.3; fewer realisations than L=12; two bands per realisation (around E and E+pi); same J, h per realisation across epsilon from a disorder table made up front; jobs split over epsilon.
- **Cost:** ~20k core-h (measured on the laptop 2026-11-10, nev=64: L=14 ~200 s per realisation and eps for both bands -> ~7 min on BAF -> ~4.7k; L=16 extrapolated ~1 h laptop -> ~2 h on BAF -> ~15k); L=14 jobs = 1 realisation x all eps of a theta (~2.5 h), L=16 jobs = 1 realisation x ~5 eps (~10 h).
- **For:** the band gap and eta crossings (Fig. 3.23 `fig:band_gap_median`, Fig. 3.20 `fig:psd_scalar`) and the time crystal boundaries of Fig. 3.24 `fig:phase_portraits_gap`, against Fig. 2.4 `fig:lsr_phase_portrait`. This decides whether the DTC boundary is distinct from the MBL one (the backbone of the thesis and of the flagship paper).
- **Notes:** the longest share of the cluster; only quantities the thesis plots.

### V1.3: rerun the five old grids that are a different circuit
- **Status:** PROPOSED 2026-11-12, waits for Tony
- **Why:** the old 0.025 grids of L=10 theta 0.05, 0.15, 0.2 and L=12 theta 0, 0.05 (rawdata/eigen, 2024) were made with the older brickwall, in which the detuning noise delh is drawn once and reused on every gate; everything else (2026-09, V1, V1.2) uses circuit_dtc with delh per gate (stored [ENV] Code checked). They are left out of the merged analysis (merged.jl OLD_MODEL), so L=12 theta 0 is only the 18 window points and L=12 theta 0.05 has nothing.
- **Parameters:** mbl_dtc set grid (eps/grid.txt, 41 points), same code as V1.2: L=12 theta 0 jobs 1-2000 (matching its windows), theta 0.05 jobs 1-1000; L=10 theta 0.05, 0.15, 0.2 at 10 realisations per job, 400 jobs each (4000 realisations, as L=10 elsewhere). Needs a new set name (e.g. grid) only where it does not collide: L=12 theta 0 and 0.05 have no grid files yet; L=10 has none.
- **Cost:** L=12 ~8.8 h per job (41 points, measured on the V1.2 grid): theta 0 ~17.6k, theta 0.05 ~8.8k core-h; L=10 ~1.4 h per job of 10: ~1.7k. Total ~28k core-h (theta 0 at L=12 alone 17.6k; 1000 jobs there would halve it).
- **For:** Fig. 2.4, 3.5, 3.25 and every theta=0 scan of Ch. 3 at L=12, and the FSS of the new section at theta 0 (today on the windows only).
