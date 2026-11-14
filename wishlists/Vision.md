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
- **Status:** approved 2026-11-10 (Tony): old theta 0, 0.1, 0.2, 0.3 first, ~20 eps per theta, L=14 x1000, L=16 x200 realisations (raised by Tony 2026-11-12: POLFED gives 128 levels per realisation, safe rather than sorry, a top-up later would wait for the priority to recover). Scripted in BAF/mbl_dtc_polfed (polfed_bands.jl, two bands at phi=0, pi, nev=64; checked against full diagonalisation at L=8, exact to 1e-14, inner half of each band complete). READY 2026-11-12 (6e6687f): tarball julia-1.11.1-polfed (Tony), eps lists 10-12 per theta, 8400 jobs; local test in the tarball's exact Manifest ok, dry run ok. BAF test skipped (Tony: keep the cluster share), queue V2 then V1.3
- **Where:** cluster (BAF)
- **Parameters:** L=14, 16 (L=18 postponed), one non-uniform epsilon list per theta, theta 0-0.3; fewer realisations than L=12; two bands per realisation (around E and E+pi); same J, h per realisation across epsilon from a disorder table made up front; jobs split over epsilon.
- **Cost:** ~40k core-h at the raised counts with ~20 eps per theta, ~24k with ~12 eps concentrated at the crossings (was ~20k at x500/x100; measured on the laptop 2026-11-10, nev=64: L=14 ~200 s per realisation and eps for both bands -> ~7 min on BAF -> ~4.7k; L=16 extrapolated ~1 h laptop -> ~2 h on BAF -> ~15k); L=14 jobs = 1 realisation x all eps of a theta (~2.5 h), L=16 jobs = 1 realisation x ~5 eps (~10 h).
- **For:** the band gap and eta crossings (Fig. 3.23 `fig:band_gap_median`, Fig. 3.20 `fig:psd_scalar`) and the time crystal boundaries of Fig. 3.24 `fig:phase_portraits_gap`, against Fig. 2.4 `fig:lsr_phase_portrait`. This decides whether the DTC boundary is distinct from the MBL one (the backbone of the thesis and of the flagship paper).
- **Notes:** the longest share of the cluster; only quantities the thesis plots.

### V1.3: rerun of the old-circuit grids at the four theta, plus the scrambling reference
- **Status:** APPROVED 2026-11-12 (Tony), scripted, mbl_dtc/jobs_v1.3.txt (4816 jobs), local test passed (julia 1.9.4, set testscr); BAF test line testscr in jobs_test.txt. Goes AFTER V2.
- **Why:** the full scan of 2026-11-12 (every file, [ENV] Code) found the older brickwall (same detuning noise on every gate) in eigen L10 theta 0.05, 0.15, 0.2 and L12 theta 0, 0.05 (all files), and in eigen_h0 L8, L10, L12 theta 0 and L12 theta 0.3 (h=0 not rerun now, kept). Everything else is circuit_dtc. Only the four theta of the thesis (0, 0.1, 0.2, 0.3) are rerun; 0.05, 0.15, 0.25, 0.35 come back with the finer theta scan. The five old random-h sets are deleted from E:.
- **Parameters:** set grid (eps/grid.txt, 41 points): L=12 theta 0, 2000 jobs x 1; L=10 theta 0.2, 400 jobs x 10. Scrambling reference theta=1, g=0.5 (set scr, eps/scr.txt): L=8 4000 (16 x 250), L=10 4000 (400 x 10), L=12 2000.
- **Cost:** ~17.6k (L=12 theta 0) + ~0.6k (L=10 theta 0.2) + ~0.5k (scr) core-h, ~19k total.
- **For:** the FSS at theta 0 (Sec. fss), every theta=0 scan of Ch. 3 at L=12, Fig. 2.4/3.5/3.25 at theta 0.2, and the scrambling examples of Ch. 3 (in place of theta=0.35, g=0.5).
- **V2 too:** theta=1, g=0.5 at L=14 x1000 and L=16 x200 (mbl_dtc_polfed/jobs.txt, eps lists ready), 1200 jobs, ~0.6k core-h.

### V2.1: L = 14, 16 on a fine epsilon grid around the crossings, enough realisations for smooth curves (kind 2)
- **Status:** APPROVED 2026-11-14 in full (Tony), L=14 and L=16 as two job lists (jobs_v2.1_L14.txt, jobs_v2.1_L16.txt, set v21, realisations after V2); with the undersampled Lazarides lists (polfed_bands.jl, Tony 2026-11-14), checked at L=8 against full diagonalisation; one-job BAF test jobs_test_v21.txt first. Planned 2026-11-14 (Tony: "that's the plan for V2.1 anyway... go for enough realisations already so that the curves are smooth enough... it's ok if it takes long to run"). Counts and grid below are my proposal, not yet approved; same script and tarball as V2 (BAF/mbl_dtc_polfed, polfed_bands.jl, nev=64), new eps lists only.
- **Why:** with V2 the FSS of <r> at L=8..16 does not resolve mu better than L<=12 (theta=0: PM 0.812(53) vs DTC 1.005(131), 1.4 sigma, as before; theta=0.1 no separation): only ~4-5 V2 points per +-0.1 window, L=16 at 200 realisations (per-realisation sd of the band <r> 0.03-0.08 -> SE 0.003), and the crossings of L=14,16 sit at the edge of the windows centred on the L<=12 crossings.
- **Where:** cluster (BAF)
- **Parameters:** theta 0 and 0.1 (where L=14,16 still cross on both sides); a uniform grid of dg=0.01 over both crossover regions: theta=0 PM g=0.14-0.28, DTC g=0.76-0.92; theta=0.1 PM g=0.05-0.22, DTC g=0.82-0.97 (about 15 points per side, ~30 eps per theta). Realisations: L=14 x2000, L=16 x1000 (SE of <r> per point ~0.001 and ~0.0015, as L=12 now). Same seeds as V2 (realisation r = the same circuit), so V2 realisations at shared points can be merged.
- **Cost (measured on V2):** L=14 ~0.73 h per (realisation, eps) for both bands, L=16 ~2.3 h (median band 4000 s, max 12000 s). L=14: 2000 x 30 x 2 x 0.73 ~ 88k core-h; L=16: 1000 x 30 x 2 x 2.3 ~ 138k core-h; ~225k core-h in all (V2 was ~65k). Cheaper variants: dg=0.015 (~20 eps per theta) ~150k; or L=16 x600 ~170k. Jobs: L=14 one realisation x 15 eps (~11 h, two jobs per realisation and theta), L=16 one (realisation, eps) per job (2-7 h), memory 2/3 GB.
- **For:** the FSS of both crossovers with L=8..16 (Fig. lsr_gc_mu, the mu comparison PM vs DTC = Tony's central message), smooth L=14,16 curves in all Ch3 scans, the undersampled PSD (eta) and band gap crossings, and F.R.I.D.A.Y.'s undersampled SFF.
- **Notes:** run after V1.3; long, so maybe in two halves (theta=0 first).

### V2.2: rerun of V2 for more realisations, now with the Lazarides lists (kind 2)
- **Status:** planned 2026-11-14 (Tony): V2's points again with more realisations; polfed_bands.jl now also stores the undersampled Lazarides lists, so these come from V2.2 only. Everything else uses V2 + V2.2 together; the Lazarides statistics only V2.2's realisations ("which I can live with"). Counts and cost to be set; realisations after V2's (and V2.1's) so they are new circuits.
- **Where:** cluster (BAF)
- **For:** the Lazarides diagrams, histograms, kurtosis and KL of Ch3 at L=14, 16 (undersampled), and more statistics at the V2 points of every Ch3 quantity.
