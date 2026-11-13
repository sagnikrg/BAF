# K.A.R.E.N. (Chapter 5): simulation wishlist

Chapter 5: entanglement structures of the DTC eigenstates (information lattice) and the Gullans-Huse type reference qubit probe ("web shooting"). Code: info_lattice_fisher_information/code/gullans_huse (gh_probe.jl, test.gh_probe.jl, script.gh_probe_small.jl), circuit_dtc from the header. Raw data on E:/info_lattice_fisher_information/rawdata/<run>/L*/ (moved 2026-11-13, git tracking rule); plotdata and logs in info_lattice_fisher_information/data/. Timings, single thread: one realisation of 100 periods with all mutual informations, 0.24 s (L=8), 1.4 s (L=10).

Existing info lattice data (mbl_dtc/code/info_lattice) is L=8 only, mostly at theta=0; eigenvectors are not stored on E:, so new info lattice figures need new runs (on hold, Tony 2026-11-08: scans and plots for the info lattice come later).

### K1: small test of the reference qubit probe
- **Status:** final (E:/info_lattice_fisher_information/rawdata/gh_probe/L8, L10), 2026-11-09; L=8 104 s, L=10 635 s
- **Where:** laptop (job queue, 1 core)
- **Parameters:** L=8 and 10, x0=L/2, T=100 periods, 50 realisations per point; the three examples (PM theta=0 g=0.05, scrambling theta=0.35 g=0.5, DTC theta=0 g=0.95) and theta=0, g=0.1:0.1:0.9
- **Cost:** about 0.3 core-h, one job
- **For:** no figure yet; first look at whether I(R:x0), I(R:block) and the R-site correlators tell DTC, scrambling and MBL apart, before Tony's workarounds for the single reference qubit
- **Notes:** quick

### K2: free Majorana scan (windings, mode splittings)
- **Status:** final (data/majorana/free_scan.hdf5), 2026-11-12
- **Where:** laptop, 1 core
- **Parameters:** free chain, L=32, 64, 128, g scan
- **For:** Fig. fig:majorana_invariants (Ch5 Sec. majorana_pi_mode)

### K3: edge autocorrelator ED (A_pi, A_0)
- **Status:** final (data/majorana/edge_ed.hdf5), 2026-11-13; L10/h1/theta0.3 redone after the Windows restart, 5 threads, 20 min
- **Where:** laptop
- **Parameters:** L=8 (100 realisations), L=10 (16); h=0 and h uniform in [0,2pi); theta 0, 0.1, 0.2, 0.3; 10 g values; t=200..400
- **For:** Fig. fig:majorana_edge

### K4: canonical driven Kitaev chain survey (converse)
- **Status:** final, 2026-11-13; raw on E:/info_lattice_fisher_information/rawdata/majorana/L*/canonical.hdf5, plotdata data/majorana/canonical_plot.hdf5 (post.canonical.jl)
- **Where:** laptop, 1 core; L=12 about 14 min per point
- **Parameters:** 8 points (0SG, PM, piSG weak/strong, 0piPM, piSG w!=Delta, cluster weak/strong); L=8, 10, 12 with 200, 30, 5 realisations
- **For:** Fig. fig:majorana_converse (Sec. majorana_converse)

### K5: A_pi at L=14, 16 with fields and interaction (idea, not approved)
- **Status:** planned, waiting for Tony (Ch5 outlook: does the dressed pi mode survive at larger L)
- **Where:** BAF; before any submission, rewrite BAF/Info_lattice_fisher_information (still the 2024 one jdl per job scripts) into the batch layout of HOWTO_cluster.md (template mbl_dtc_polfed/cat_batch.sh), interactive test, dry run
- **Parameters:** g = 0.9, 0.95, 1; theta 0-0.3; h=0 and random; L=12, 14, 16
- **Cost:** not estimated yet

### K6: reference qubit probe at the example and Kitaev points
- **Status:** final (data/majorana/gh_edge.hdf5), 2026-11-13, about 7 min on 1 core
- **Parameters:** L=8 (50), L=10 (30) realisations; R on edge and bulk qubit; t=100..200
- **For:** Fig. fig:gh_probe_classes (Ch5 Sec. gh_probe)

### K7: reference qubit probe g scans
- **Status:** final (data/majorana/gh_scan.hdf5), 2026-11-13; 6 threads at the office, L=10 theta=0.4 single thread at home (27 min)
- **Parameters:** g = 0:0.05:1, theta = 0..0.4, fields on; L=8 (50), L=10 (30)
- **For:** Fig. fig:gh_probe_scan; L=12+ (time evolution only) is the natural next run for the probe chapter
