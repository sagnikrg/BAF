# E.D.I.T.H. (Chapter 4, Vertical Phases): simulation wishlist

Figure numbers from the full Bonn draft of 2026-11-08 (14:35). Random h, the noisy circuit of the thesis, vertical_phases code (`code/run_vertical.jl`, copy in `BAF/vertical_phases/code`). One realisation = U_0 (theta=0, same disorder) plus U_F at every theta of the scan, all dense; default theta grid 0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.5, 0.75, 1 (13 values). Jobs are generated from `BAF/vertical_phases/jobs.txt`. Timings are from one laptop core: L=10 ~25 s, L=12 ~26 min per realisation.

### E1: night run, g = 0.25, 0.5, 0.75, theta 0-1 (8 values)
- **Status:** final (`vertical_phases/data/vertical_unitary/`), 2026-11-08
- **Where:** laptop (job queue, vp_night)
- **Parameters:** theta 0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1; g=0.5: L=8, 10, 12 with 200, 100, 30; g=0.25, 0.75: L=8, 10 with 200, 60
- **Cost:** 9.3 h on one core
- **For:** all Chapter 4 figures as they are in the draft now (Fig. 4.1-4.4); replaced by E2/E3 once those are back
- **Notes:** old file layout, same seeds as E2/E3, never mix the two in one data folder

### E2a: g = 0.5, fine theta grid, L = 8 and 10 (QUICK)
- **Status:** final (raw `E:/vertical_phases`, reduced `vertical_phases/plotdata/vertical_cluster.hdf5`), 2026-11-10; 13 of the L=10 jobs missing (L=10: 600/1400/700 of 1000/2000/1000 at g=0.25/0.5/0.75), enough for the figures
- **Where:** cluster (BAF)
- **Parameters:** g=0.5, 13 thetas; L=8 5000 realisations (jobs of 1000), L=10 2000 (jobs of 100)
- **Cost:** ~17 core-h; 25 jobs, longest ~45 min
- **For:** Fig. 4.1 `fig:vertical_V_G`, Fig. 4.2 `fig:vertical_G_res` (panel b.), Fig. 4.4 `fig:vertical_measures`, and the new figure of the peak at delta=0 on log bins (`vertical_peak.png`, not in the draft yet). Settles the shape of the rise between theta=0 and 0.35 and the bands at small theta, where the mean of |V|^2 is dominated by rare near resonances
- **Notes:** quick, can go any time

### E2b: g = 0.5, fine theta grid, L = 12
- **Status:** final (raw `E:/vertical_phases`, reduced `vertical_phases/plotdata/vertical_cluster.hdf5`), 2026-11-10; 13 of the L=10 jobs missing (L=10: 600/1400/700 of 1000/2000/1000 at g=0.25/0.5/0.75), enough for the figures
- **Where:** cluster (BAF)
- **Parameters:** g=0.5, 13 thetas, L=12, 1000 realisations (jobs of 2)
- **Cost:** ~430 core-h; 500 jobs of ~1 h, 6 GB each
- **For:** same figures as E2a, plus Fig. 4.3 `fig:vertical_pdelta` (L=12 only); the "bigger L scrambles more" point and the absence of crossings in L rest on these curves
- **Notes:** many short jobs, can fill gaps between long ones

### E3a: g = 0.25 and 0.75, L = 8 and 10 (QUICK)
- **Status:** final (raw `E:/vertical_phases`, reduced `vertical_phases/plotdata/vertical_cluster.hdf5`), 2026-11-10; 13 of the L=10 jobs missing (L=10: 600/1400/700 of 1000/2000/1000 at g=0.25/0.5/0.75), enough for the figures
- **Where:** cluster (BAF)
- **Parameters:** g=0.25 and 0.75, 13 thetas; L=8 2000 each (jobs of 1000), L=10 1000 each (jobs of 100)
- **Cost:** ~15 core-h; 24 jobs, longest ~45 min
- **For:** Fig. 4.2 `fig:vertical_G_res` (panels a. and c.), and the sentence that the other g behave the same (optionally one extra figure of the four measures at g=0.25, 0.75)
- **Notes:** quick

### E3b: g = 0.25 and 0.75, L = 12
- **Status:** final (raw `E:/vertical_phases`, reduced `vertical_phases/plotdata/vertical_cluster.hdf5`), 2026-11-10; 13 of the L=10 jobs missing (L=10: 600/1400/700 of 1000/2000/1000 at g=0.25/0.5/0.75), enough for the figures
- **Where:** cluster (BAF)
- **Parameters:** g=0.25 and 0.75, 13 thetas, L=12, 500 realisations each (jobs of 2)
- **Cost:** ~430 core-h; 500 jobs of ~1 h, 6 GB each
- **For:** Fig. 4.2 `fig:vertical_G_res` panels a. and c. get a third size
- **Notes:** lower priority than E2b

### E4: L = 14 at g = 0.5
- **Status:** planned (held for later by Tony, 2026-11-08); commented out in jobs.txt
- **Where:** cluster (BAF)
- **Parameters:** L=14, g=0.5, theta 0, 0.05, 0.1, 0.2, 0.35, 0.5, 1; 10 realisations, one per job, 8 threads
- **Cost:** estimate ~5 h per job on 8 threads (~40 core-h each, ~400 core-h total), 48 GB per job; 10 jobs
- **For:** one more size in Fig. 4.2-4.4 for the finite size part of the key idea
- **Notes:** memory heavy (~4 GB per dense matrix); after the L=12 runs

### E5: the 13 missing L = 10 jobs of E2/E3 (QUICK)
- **Status:** planned, 2026-11-10
- **Where:** cluster (BAF)
- **Parameters:** the L=10 blocks of 100 not in E:/vertical_phases (g=0.25: 4, g=0.5: 6, g=0.75: 3); same seeds, so they complete the planned 1000/2000/1000
- **Cost:** ~9 core-h; 13 jobs of ~40 min
- **For:** Fig. 4.2-4.5 (`fig:vertical_G_res`, `fig:vertical_pdelta`, `fig:vertical_peak`, `fig:vertical_measures`), L=10 bands only
- **Notes:** quick, needs a jobs list with just those ranges

### E6: check of the G_res ceiling with a Haar random drive (QUICK)
- **Status:** planned, 2026-11-10
- **Where:** laptop (job queue)
- **Parameters:** U_0 of the real circuit at g=0.5, V replaced by a Haar random unitary; L=8, 10, 12 with 200, 100, 20 realisations; same estimator (window |dphi|<2pi/96)
- **Cost:** ~1 core-h
- **For:** Fig. 4.3 `fig:vertical_G_res`: whether the 10% overshoot of L=12 above 2^L/(2pi)^2 at theta=0.5-0.75 is the near degenerate spike at dphi->0 (Fig. 4.2 `fig:vertical_V_G`) or the estimator
- **Notes:** settles one sentence in the text

### E7: coarse (g, theta) map at L = 8, 10 (optional)
- **Status:** planned, 2026-11-10 (Outlook material; Tony: weight competition stays an idea in the thesis)
- **Where:** cluster (BAF)
- **Parameters:** g=0.05..0.95 in steps of 0.1, the 13 thetas; L=8 500, L=10 200 realisations
- **Cost:** ~16 core-h; ~30 jobs, longest ~40 min
- **For:** the Outlook `sec:vertical_outlook` (largest overlap and S_half over the (g, theta) plane); not a thesis figure unless Tony wants one
- **Notes:** quick
