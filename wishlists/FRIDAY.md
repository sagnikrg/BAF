# F.R.I.D.A.Y. (Chapter 6, Spectral Form Factor): simulation wishlist

Figure numbers from the full draft of 2026-11-08 (main.aux 14:15). All figures of Chapter 6 (6.1-6.9) are complete at L = 8, 10, 12 from the stored eigenvalues on E:. **Nothing new is needed for the chapter as it stands.** The two planned entries below are cheap, and they only matter if L = 14, 16 or the envelope scaling go into the text.

### F1: envelope on a t/N grid, L = 8, 10, 12, theta = 0-0.3
- **Status:** final (`mbl_dtc/code/eigen_signatures/boot/sff_env_L{8,10,12}.hdf5`), 2026-11-08
- **Where:** laptop (jobs sff_env_tau_L8, sff_env_tau_L10, sff_env_tau_L12)
- **Parameters:** g = 0.025:0.025:0.975, theta = 0, 0.1, 0.2, 0.3, tau = t/N = 2^(0:0.5:6), window of 200 steps. Realisations at theta = 0: 8000 / 4000 / 4066; at 0.1: 4000 / 2000 / 500; at 0.2: 4000 / 3920 / 500; at 0.3: 4000 / 2000 / 500 (L = 8 / 10 / 12). Also the histograms of ln s, where s is the pi-pair splitting in units of the level spacing.
- **Cost:** about 3.5 core-h in total (45 min, 64 min and 127 min, single threaded, E: reader)
- **For:** Fig. 6.4 `fig:sff_envelope` and Fig. 6.5 `fig:sff_envelope_theta` (the tau = 4 slice). The tau grid supports the scaling analysis in `mbl_dtc/code/eigen_signatures/notes/envelope_scaling.md`, which is not in the text yet.

### F2: dense splitting check of the clean circuit, theta = 0, L = 6, 8, 10
- **Status:** final (`mbl_dtc/code/eigen_signatures/notes/splitting_check.log`), 2026-11-07
- **Where:** laptop (job split_check)
- **Parameters:** clean circuit (no tuning noise), theta = 0, g = 0.70-0.95, 400 / 200 / 60 realisations at L = 6 / 8 / 10
- **Cost:** under 0.1 core-h
- **For:** no figure. It checks var(ln s) proportional to L and the crossing g_c of about 0.85 in the scaling note.

### F3: dense splitting check of the clean circuit at L = 12 (quick)
- **Status:** planned, 2026-11-08
- **Where:** laptop (job queue)
- **Parameters:** clean circuit, theta = 0, g = 0.75-0.975 in steps of 0.025 (10 points), about 100 realisations per point, L = 12
- **Cost:** about 10-20 core-h, one single-threaded job, no E: needed
- **For:** no figure yet. It fixes lambda(epsilon) and g_c without the noise floor, which the raw data at L = 12 has (mean ln s of about -11.7 against about -15 expected clean at g = 0.95). Only needed if the t/N scaling result (the stability criterion lambda(epsilon) < 0) goes into Chapter 6.
- **Notes:** quick, and the only real run on this list.

### F4: L = 14, 16 envelope, riding on Vision's V2 (quick, zero extra runs)
- **Status:** planned, 2026-11-08
- **Where:** cluster (BAF), inside V2
- **Parameters:** whatever V2 runs. The one request is that V2 keeps, per realisation and point, both antipodal bands (it already plans to), so that the distance of each level to its partner at phi + pi can be read off.
- **Cost:** zero extra core-h. Reading the splittings afterwards is minutes on the laptop.
- **For:** Fig. 6.4 `fig:sff_envelope`, the extension to L = 14, 16. At late times E = 2<cos 2 pi tau s> holds to 0.01 at tau >= 8 (checked at L = 8, 10, 12), so the splittings alone give the envelope without the full spectrum. At tau = 4 the agreement is only 0.02-0.08, so the figure would move to tau = 8. This settles whether E at fixed t/N keeps growing with L in the time crystal (it does from 8 to 12 at g >= 0.85).
- **Notes:** depends on V2. Nothing to submit separately.

### Not needed
- The old pinned rerun of the SFF envelope at L = 12, theta = 0.05 (files that failed in the drive drop of 2026-11-06): theta = 0.05 is not in the thesis anymore, so drop it.
