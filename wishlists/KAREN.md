# K.A.R.E.N. (Chapter 5): simulation wishlist

Chapter 5: entanglement structures of the DTC eigenstates (information lattice) and the Gullans-Huse type reference qubit probe ("web shooting"). Code: info_lattice_fisher_information/code/gullans_huse (gh_probe.jl, test.gh_probe.jl, script.gh_probe_small.jl), circuit_dtc from the header. Data: info_lattice_fisher_information/data/gh_probe, never E:. Timings, single thread: one realisation of 100 periods with all mutual informations, 0.24 s (L=8), 1.4 s (L=10).

Existing info lattice data (mbl_dtc/code/info_lattice) is L=8 only, mostly at theta=0; eigenvectors are not stored on E:, so new info lattice figures need new runs (on hold, Tony 2026-11-08: scans and plots for the info lattice come later).

### K1: small test of the reference qubit probe
- **Status:** final (info_lattice_fisher_information/data/gh_probe/gh_probe_small_L8.hdf5, _L10.hdf5), 2026-11-09; L=8 104 s, L=10 635 s
- **Where:** laptop (job queue, 1 core)
- **Parameters:** L=8 and 10, x0=L/2, T=100 periods, 50 realisations per point; the three examples (PM theta=0 g=0.05, scrambling theta=0.35 g=0.5, DTC theta=0 g=0.95) and theta=0, g=0.1:0.1:0.9
- **Cost:** about 0.3 core-h, one job
- **For:** no figure yet; first look at whether I(R:x0), I(R:block) and the R-site correlators tell DTC, scrambling and MBL apart, before Tony's workarounds for the single reference qubit
- **Notes:** quick
