# J.A.R.V.I.S. (Chapter 7): simulation wishlist

Figure numbers from the full draft of 2026-11-08. Random h, the noisy circuit of the thesis. Code: sff/code/symmetry_breaking (sb.jl exact diagonalisation of U_F, open and periodic chains; tm.jl short time transfer matrix). Costs are scaled from laptop timings, single thread: diagonalisation without vectors 2.7 s (L=10) and 106 s (L=12), with vectors and order parameters about 165 s (L=12), L=14 extrapolated at about 1.3 h. Output is small (K(t) arrays), goes to sff/rawdata/sb, never to E:.

Quick ones: J1, J2 (laptop); J3a, J5, J6, J7 (cluster, each at most about 200 core-h).

### J1: restart of sb_ch7_night (run.sb.jl), the last L=12 point
- **Status:** approved (night of 2026-11-07/08), interrupted by the pause at 16:51 on 2026-11-08; restarts when Tony is home. Skips finished tasks.
- **Where:** laptop (job queue, 1 core)
- **Parameters:** L=12, theta=0, g=0.95, open and periodic chain, 40 realisations (K only)
- **Cost:** about 2 core-h, one job
- **For:** the L=12 point of the planned lambda(0,t) and count figure (`fig:gc_lambda_count`, not in the draft yet), see J3
- **Notes:** quick; on its own too few realisations for the figure, J3 replaces it

### J2: transfer matrix at t=3 (tm_run.jl 3 200)
- **Status:** approved (planned after J1), 2026-11-08
- **Where:** laptop (job queue, 1 core)
- **Parameters:** t=3, theta 0 and 0.35, g 0:0.05:1 (42 points), 200 noise samples; exact averages over J and h
- **Cost:** about 1 core-h, one job
- **For:** Fig. 7.1 `fig:gc_tm_spectra`, the empty right column. Settles the first time where the time crystal and the scrambling phase differ: at odd t the time crystal has no t/2 sector, so a second leading eigenvalue there means scrambling, none means DTC.
- **Notes:** quick

### J3: lambda(0,t) and the count of leading eigenvalues, L = 8, 10, 12
- **Status:** planned, 2026-11-08
- **Where:** cluster (BAF)
- **Parameters:** theta=0 at g=0.05:0.1:0.95, plus theta=0.35, g=0.5 (12 points); open and periodic chain, eigenvalues only; realisations 20000 (L=8), 10000 (L=10), 4000 (L=12)
- **Cost:** J3a (L=8, 10): about 230 core-h, 500 realisations per job, about 45 min each, about 300 jobs. J3b (L=12): about 2800 core-h, 100 realisations per job, about 6 h each, about 480 jobs.
- **For:** a new figure `fig:gc_lambda_count` (the sharpest test of Garratt and Chalker): lambda(0,t) from the L scaling of K_open, and K_periodic/lambda^L counting the leading eigenvalues, which should be 1 (paramagnet), 2 at even t (time crystal) and t (scrambling). With 300 realisations at L=10 it is pure noise today.
- **Notes:** J3a is quick and already shows whether J3b is needed at full size

### J4: the same at L = 14
- **Status:** planned, 2026-11-08
- **Where:** cluster (BAF)
- **Parameters:** L=14, the three examples (theta 0 g 0.05; theta 0.35 g 0.5; theta 0 g 0.95), open and periodic, 200 realisations
- **Cost:** about 1600 core-h, one realisation per job, about 2.6 h each, 600 jobs, about 10 GB RAM per job
- **For:** `fig:gc_lambda_count`, a fourth size for the fit
- **Notes:** long; only if J3 is not convincing with three sizes

### J5: response to fields of period b, L = 10, 12
- **Status:** planned, 2026-11-08
- **Where:** cluster (BAF)
- **Parameters:** the three examples, b = 2, 3, 4, field on one site; 2000 realisations (L=10), 500 (L=12)
- **Cost:** about 200 core-h, about 100 jobs of about 2 h
- **For:** Fig. 7.4 `fig:gc_kb`. L=8 already shows K_b/K = 1/gcd(b,2) in the time crystal; L=10 (150 realisations today) is too noisy for bands. Settles that the parity response holds with size, and whether the scrambling point reaches 1/b.
- **Notes:** quick

### J6: order parameters at L = 12, and an L = 10 g scan
- **Status:** planned, 2026-11-08
- **Where:** cluster (BAF)
- **Parameters:** L=12 at the three examples, 1000 realisations; L=10 at theta=0, g=0.05:0.1:0.95, 1000 realisations; eigenvectors, O2 for nu=1 and nu=t/2
- **Cost:** about 170 core-h, about 100 jobs of 1-2 h
- **For:** Fig. 7.2 `fig:gc_o2_profiles` (third size) and Fig. 7.3 `fig:gc_o2_time` (L=10 next to L=8). Settles the long range order of O2^(t/2) in the time crystal against the decay in the paramagnet with size.
- **Notes:** quick

### J7: a cleaner ergodic reference point
- **Status:** planned, 2026-11-08
- **Where:** cluster (BAF)
- **Parameters:** g=0.5, theta = 0.5 and 0.75, L=10 (1000 realisations) and L=12 (300), order parameters and fields of period b
- **Cost:** about 100 core-h, about 60 jobs
- **For:** the scrambling column of Figs. 7.2, 7.4 (`fig:gc_o2_profiles`, `fig:gc_kb`). At theta=0.35 the approach to the ergodic values is slow (O2^(t/2)/t^2 = 0.22 at L=8, t=64; K_b/K still drifting to 1/b). Settles whether that is the point or the sizes.
- **Notes:** quick. Theta=1 stays on the after-thesis list.

### J8: exact transfer matrix at t = 4, 5 (no finite size)
- **Status:** planned, 2026-11-08; needs new code (matrix free, J and h averages as exact discrete Fourier sums)
- **Where:** laptop for t=4, cluster for t=5
- **Parameters:** theta 0 and 0.35, g grid 0:0.05:1, about 200 noise samples
- **Cost:** not known until the code exists
- **For:** Fig. 7.1 `fig:gc_tm_spectra`, more columns. Leading eigenvalues per sector of the infinite chain at t=4 and the odd t=5: the mod two of the time crystal in its cleanest form.
- **Notes:** the strongest figure of the chapter if it works
